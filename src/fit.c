#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linesearch.h"
#include "fit.h"


// sets model and loglik
// advances frame to end of messages
// assumes coefs has already been set
static void evaluate_loglik(struct recv_fit *fit)
{
	frame_clear(&fit->frame);
	
	// update coefs
	model_set_recv_coefs(&fit->model, &fit->coefs);

	// update loglik
	recv_loglik_clear(&fit->loglik);
	recv_loglik_add_all(&fit->loglik, &fit->frame, fit->msgs);
}

// sets scale
// assumes evaluate has already been called
static void compute_scale(struct recv_fit *fit)
{
	const double var_tol = 1e-8;	
	ssize_t i, dim = vector_dim(&fit->coefs);
	
	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik);
	const struct matrix *imat = &info->imat;
	
	/* get the empirical variances of the covariates */
	assert(vector_dim(&fit->scale) == matrix_diag_dim(imat, 0));
	matrix_get_diag(imat, 0, vector_to_ptr(&fit->scale));

	/* define the scale to be 1 when the variance is 0 */
	for (i = 0; i < dim; i++) {
		double *s = vector_item_ptr(&fit->scale, i);
		if (*s < var_tol)
			*s = 1;
	}
	
	/* otherwise, take the square root of the variance */
	vector_sqrt(&fit->scale);
	vector_fill(&fit->scale, 1.0);
}

// sets ce_t, be, and ne
// reinits kkt and resid, duals, and ldlfac
// assumes scale is initalized and evaluate has already been called
static void compute_constraints(struct recv_fit *fit)
{
	const double eig_tol = 1e-5;
	ssize_t i, dim = vector_dim(&fit->coefs);
	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik);
	const struct matrix *imat = &info->imat;
	
	/* compute the eigendecomp of the empirical correlation matrix;
	 * kkt matrix and resid vector are un-initialized/free, so we
	 * can use them for temporary storage.
	 */
	struct matrix evec = matrix_slice(&fit->kkt, 0, 0, dim, dim);
	struct vector eval = vector_slice(&fit->resid, 0, dim);

	matrix_assign_copy(&evec, imat);
	matrix_div_rows(&evec, &fit->scale);
	matrix_div_cols(&evec, &fit->scale);

	/* compute its eigendecomposition */
	assert(symeig_dim(&fit->eig) == dim);
	assert(symeig_job(&fit->eig) == EIG_VEC);
	bool ok = symeig_factor(&fit->eig, UPLO_LOWER, &evec, &eval);
	assert(ok);
	
	/* compute its rank */
	for (i = 0; i < dim; i++) {
		double l = vector_item(&eval, i);
		if (l > eig_tol)
			break;
	}
	// rank = dim - i;
	
	/* define the constraints */
	struct matrix u = matrix_slice_cols(&evec, 0, i);
	matrix_reinit(&fit->ce_t, dim, i);
	matrix_assign_copy(&fit->ce_t, &u);
	vector_reinit(&fit->be, i);
	vector_fill(&fit->be, 0);
	fit->ne = i;
	
	/* update dims of kkt matrix, residuals and dual parameters */
	matrix_reinit(&fit->kkt, dim + i, dim + i);
	vector_reinit(&fit->resid, dim + i);
	vector_reinit(&fit->params, dim + i);
	
	/* update ldlfac dimension */
	ldlfac_reinit(&fit->ldl, dim + i);
}

// sets resid
// assumes evaluate and preprocess have already been called
static void compute_resid(struct recv_fit *fit)
{
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = fit->ne;

	assert(vector_dim(&fit->be) == ne);
	
	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik);
	const struct vector *score = &info->score;
	
	/* r1 is the dual residual: grad(f) + ce' * nu,
	 * where grad(f) = s^{-1} * grad(nll) + lambda * x
	 *               = -[s^{-1} * score - lambda * x]
	 */
	struct vector r1 = vector_slice(&fit->resid, 0, dim);
	struct vector primals = vector_slice(&fit->params, 0, dim);
	struct vector duals = vector_slice(&fit->params, dim, ne);
	vector_assign_copy(&r1, score);
	vector_div(&r1, &fit->scale);
	vector_axpy(-fit->penalty, &primals, &r1);
	matrix_mul(1.0, TRANS_NOTRANS, &fit->ce_t, &duals, -1.0, &r1);
	
	/* r2 is the primal residual: ce * x - be */
	struct vector r2 = vector_slice(&fit->resid, dim, ne);
	vector_assign_copy(&r2, &fit->be);
	matrix_mul(1.0, TRANS_TRANS, &fit->ce_t, &primals, -1.0, &r2);
}

// sets kkt
// assumes evaluate and preproeces have already been called
static void compute_kkt(struct recv_fit *fit)
{
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = fit->ne;

	assert(matrix_nrow(&fit->kkt) == dim + ne);
	assert(matrix_ncol(&fit->kkt) == dim + ne);

	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik);
	const struct matrix *imat = &info->imat;

	
	/* k11 is the (scaled) hessian */
	struct matrix k11 = matrix_slice(&fit->kkt, 0, 0, dim, dim);
	matrix_assign_copy(&k11, imat);
	matrix_div_rows(&k11, &fit->scale);
	matrix_div_cols(&k11, &fit->scale);
	ssize_t i;
	for (i = 0; i < dim; i++) {
		*matrix_item_ptr(&k11, i, i) += fit->penalty;
	}
	
	/* k12 is the transpose of the equality constraint matrix */
	struct matrix k12 = matrix_slice(&fit->kkt, 0, dim, dim, ne);
	matrix_assign_copy(&k12, &fit->ce_t);
	
	/* k22 is zero */
	struct matrix k22 = matrix_slice(&fit->kkt, dim, dim, ne, ne);
	matrix_fill(&k22, 0.0);
	
	/* k21 is never referenced */
#ifndef NDEBUG
	struct matrix k21 = matrix_slice(&fit->kkt, dim, 0, ne, dim);
	matrix_fill(&k21, 0);
#endif
}

// sets scale, ce, be
// reinits kkt and nresid
// destroys model and coefs 
static void preprocess(struct recv_fit *fit)
{
	/* compute the empirical covariance matrix of the covariates
	 * by using the null model (coefs = 0, all receivers equally
	 * likely
	 */
	vector_fill(&fit->coefs, 0);
	evaluate_loglik(fit);
	compute_scale(fit);
	compute_constraints(fit);
	
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = fit->ne;
	struct vector primals = vector_slice(&fit->params, 0, dim);
	struct vector duals = vector_slice(&fit->params, dim, ne);
	
	//vector_assign_copy(&fit->coefs, &fit->loglik.info.score);
	//vector_scale(&fit->coefs, 0.1/vector_max_abs(&fit->coefs));
	//vector_axpy(1.0, &fit->loglik.info.mean, &fit->coefs);
	//vector_shift(&fit->coefs, 1.0);
	//vector_log(&fit->coefs); 
	//vector_scale(&fit->coefs, 0.5); // log(sqrt(1 + obs))
	
	vector_fill(&fit->coefs, 0.0);
	
	vector_assign_copy(&primals, &fit->coefs);
	vector_mul(&primals, &fit->scale);

	vector_fill(&duals, 0.0);
	
	//vector_fill(&primals, 0.0);
	//vector_fill(&duals, 0.0);
	//vector_fill(&fit->coefs, 0.0);
	//vector_div(&fit->coefs, &fit->scale);
	evaluate_loglik(fit);
	compute_resid(fit);
	compute_kkt(fit);
}



static bool primal_dual_step(struct recv_fit *fit)
{
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = fit->ne;
	
	double f0 = vector_norm2(&fit->resid);
	double g0 = -2 * f0;

	if (g0 == 0)
		return false; /* residual is 0 */
	
	/* determine the search direction */
	struct vector search;
	vector_init_copy(&search, &fit->resid);
	struct matrix smat = matrix_make(&search, dim + ne, 1);
	ssize_t info = ldlfac_solve(&fit->ldl, UPLO_UPPER, &fit->kkt, &smat);
	assert(info == 0);
	vector_scale(&search, -1.0);
	
	struct vector grad;
	vector_init(&grad, dim + ne); /* gradients of |resid|^2 */
	struct vector cgrad = vector_slice(&grad, 0, dim);
	struct vector dgrad = vector_slice(&grad, dim, ne);

	/* save the initial primal and dual variables */
	struct vector params0;
	vector_init_copy(&params0, &fit->params);
	
	/* define the relevant subvectors */
	struct vector primals = vector_slice(&fit->params, 0, dim);
	
	/* perform a linesearch to reduce the residual norm */
	double f, g, stp;
	enum linesearch_task task;
	
	linesearch_start(&fit->ls, 1.0, f0, g0, &fit->lsctrl);
	
	printf("f0: %.22f g0: %.22f\n", f0, g0);
	
	do {
		/* take a step to get new primal and dual variables */
		stp = linesearch_step(&fit->ls);		
		vector_assign_copy(&fit->params, &params0);
		vector_axpy(stp, &search, &fit->params);

		/* compute the fit, residuals, and kkt matrix at the new point */
		vector_assign_copy(&fit->coefs, &primals);
		vector_div(&fit->coefs, &fit->scale);
		evaluate_loglik(fit);
		compute_resid(fit);
		compute_kkt(fit);
	
		struct vector cresid = vector_slice(&fit->resid, 0, dim);
		struct vector dresid = vector_slice(&fit->resid, dim, ne);
		struct matrix k11 = matrix_slice(&fit->kkt, 0, 0, dim, dim);

		/* compute the gradient of the squared residual norm */
		matrix_mul(1.0, TRANS_NOTRANS, &k11, &cresid, 0.0, &cgrad);
		matrix_mul(1.0, TRANS_NOTRANS, &fit->ce_t, &dresid, 1.0, &cgrad);
		matrix_mul(1.0, TRANS_TRANS, &fit->ce_t, &cresid, 1.0, &dgrad);
		vector_scale(&cgrad, 2.0);
		vector_scale(&dgrad, 2.0);	
		
		/* evaluate the new squared residual norm and directional derivative */
		f = vector_norm2(&fit->resid);
		g = vector_dot(&search, &grad);
		task = linesearch_advance(&fit->ls, f, g);
		
		printf("f: %.22f g: %.22f  stp: %.22f\n", f, g, stp);
	} while (task == LINESEARCH_STEP && !linesearch_sdec(&fit->ls));

	assert(linesearch_sdec(&fit->ls));
	
	vector_deinit(&params0);
	vector_deinit(&grad);	
	vector_deinit(&search);
	
	return linesearch_sdec(&fit->ls);
}


void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct vector *coefs0,
		   double penalty)
{
	assert(fit);
	assert(msgs);
	assert(design);	
	assert(messages_max_to(msgs) < design_send_count(design));
	assert(messages_max_from(msgs) < design_recv_count(design));
	assert(!coefs0 || vector_dim(coefs0) == design_recv_dim(design));
	assert(penalty >= 0.0 && isfinite(penalty));
	
	ssize_t dim = design_recv_dim(design);

	fit->design = design;
	fit->msgs = msgs;
	fit->penalty = penalty;
	
	frame_init(&fit->frame, design);		
	
	if (!coefs0) {
		vector_init(&fit->coefs, dim);
	} else {
		vector_init_copy(&fit->coefs, coefs0);
	}
	

	model_init(&fit->model, &fit->frame, &fit->coefs);
	recv_loglik_init(&fit->loglik, &fit->model);
	
	vector_init(&fit->scale, dim);
	matrix_init(&fit->ce_t, dim, 0);
	vector_init(&fit->be, 0);
	vector_init(&fit->params, dim);
	vector_init(&fit->resid, dim);
	matrix_init(&fit->kkt, dim, dim);
	symeig_init(&fit->eig, dim, EIG_VEC);
	ldlfac_init(&fit->ldl, dim);

	struct linesearch_ctrl lsctrl = LINESEARCH_CTRL0;
	fit->lsctrl = lsctrl;
	
	preprocess(fit);
}

void recv_fit_deinit(struct recv_fit *fit)
{
	assert(fit);
	
	ldlfac_deinit(&fit->ldl);
	symeig_deinit(&fit->eig);
	matrix_deinit(&fit->kkt);
	vector_deinit(&fit->resid);
	vector_deinit(&fit->params);
	vector_deinit(&fit->be);
	matrix_deinit(&fit->ce_t);	
	vector_deinit(&fit->scale);
	recv_loglik_deinit(&fit->loglik);
	model_deinit(&fit->model);
	vector_deinit(&fit->coefs);
	frame_deinit(&fit->frame);
}

bool recv_fit_step(struct recv_fit *fit)
{
	assert(fit);
	return primal_dual_step(fit);
}

bool recv_fit_converged(const struct recv_fit *fit)
{
	return false;
}

const char *recv_fit_errmsg(const struct recv_fit *fit)
{
	return NULL;
}
