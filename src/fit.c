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
	model_clear(&fit->model);
	
	// update coefs
	model_set(&fit->model, &fit->frame, &fit->coefs);

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

// sets ce_t, and be
// reinits kkt and resid, and duals
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
	matrix_scale_rows(&evec, &fit->scale);
	matrix_div_cols(&evec, &fit->scale);
	matrix_div_rows(&evec, &fit->scale);
	
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
	
	/* update dims of kkt matrix, residuals and dual parameters */
	matrix_reinit(&fit->kkt, dim + i, dim + i);
	vector_reinit(&fit->resid, dim + i);
	vector_reinit(&fit->duals, i);
}

// sets resid
// assumes evaluate and preprocess have already been called
static void compute_resid(struct recv_fit *fit)
{
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = matrix_ncol(&fit->ce_t);

	assert(vector_dim(&fit->be) == ne);
	
	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik);
	const struct vector *score = &info->score;
	
	/* r1 is the dual residual: grad(f) + ce' * nu,
	 * where grad(f) = s^{-1} * grad(nll)
	 */
	struct vector r1 = vector_slice(&fit->resid, 0, dim);
	vector_assign_copy(&r1, score);
	vector_div(&r1, &fit->scale);
	matrix_mul(1.0, TRANS_NOTRANS, &fit->ce_t, &fit->duals, -1.0, &r1);
	
	/* r2 is the primal residual: ce * x - be */
	struct vector r2 = vector_slice(&fit->resid, dim, ne);
	vector_assign_copy(&r2, &fit->be);
	matrix_mul(1.0, TRANS_TRANS, &fit->ce_t, &fit->coefs, -1.0, &r2);
}

// sets kkt
// assumes evaluate and preproeces have already been called
static void compute_kkt(struct recv_fit *fit)
{
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = matrix_ncol(&fit->ce_t);
	
	assert(matrix_nrow(&fit->kkt) == dim + ne);
	assert(matrix_ncol(&fit->kkt) == dim + ne);

	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik);
	const struct matrix *imat = &info->imat;

	/* k11 is the (scaled) hessian */
	struct matrix k11 = matrix_slice(&fit->kkt, 0, 0, dim, dim);
	matrix_assign_copy(&k11, imat);
	matrix_div_rows(&k11, &fit->scale);
	matrix_div_cols(&k11, &fit->scale);	
	
	/* k12 is the transpose of the equality constraint matrix */
	struct matrix k12 = matrix_slice(&fit->kkt, 0, dim, dim, ne);
	matrix_assign_copy(&k12, &fit->ce_t);
	
	/* k22 is zero */
	struct matrix k22 = matrix_slice(&fit->kkt, dim, dim, ne, ne);
	matrix_fill(&k22, 0.0);
	
	/* k21 is never referenced */
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
	compute_resid(fit);
	compute_kkt(fit);
}



static bool primal_dual_step(struct recv_fit *fit)
{
	ssize_t dim = model_recv_dim(&fit->model);
	ssize_t ne = vector_dim(&fit->be);
	
	double f0 = vector_norm2(&fit->resid);
	double g0 = -2 * f0;

	if (g0 == 0)
		return false; /* residual is 0 */
	
	/* determine the search direction */
	struct matrix resid_mat = matrix_make(&fit->resid, dim + ne, 1);
	ssize_t info = chol_solve(UPLO_UPPER, &fit->kkt, &resid_mat);
	assert(info == 0);
	vector_scale(&fit->resid, -1.0);
	
	struct vector csearch = vector_slice(&fit->resid, 0, dim);
	vector_mul(&csearch, &fit->scale);

	struct vector dsearch = vector_slice(&fit->resid, dim, ne);
	struct vector cgrad, dgrad; /* gradients of |resid|^2 */
	vector_init(&cgrad, dim);
	vector_init(&dgrad, ne);

	/* save the initial primal and dual variables */
	struct vector coefs0, duals0;
	vector_init_copy(&coefs0, &fit->coefs);
	vector_init_copy(&duals0, &fit->duals);
	
	/* perform a linesearch to reduce the residual norm */
	double f, g, stp;
	enum linesearch_task task;
	
	linesearch_start(&fit->ls, 1.0, f0, g0, &fit->lsctrl);
	
	printf("f0: %.22f g0: %.22f\n", f0, g0);
	
	do {
		/* take a step to get new primal and dual variables */
		stp = linesearch_step(&fit->ls);		
		vector_assign_copy(&fit->coefs, &coefs0);
		vector_assign_copy(&fit->duals, &duals0);
		vector_axpy(stp, &csearch, &fit->coefs);
		vector_axpy(stp, &dsearch, &fit->duals);

		/* compute the fit, residuals, and kkt matrix at the new point */
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
		g = vector_dot(&csearch, &cgrad) + vector_dot(&dsearch, &dgrad);
		task = linesearch_advance(&fit->ls, f, g);
		
		printf("f: %.22f g: %.22f  stp: %.22f\n", f, g, stp);
	} while (task == LINESEARCH_STEP);

	vector_deinit(&duals0);
	vector_deinit(&coefs0);
	vector_deinit(&dgrad);	
	vector_deinit(&cgrad);
	
	return task == LINESEARCH_CONV;
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
	

	model_init(&fit->model, fit->design, &fit->coefs);
	recv_loglik_init(&fit->loglik, &fit->model);
	
	vector_init(&fit->scale, dim);
	matrix_init(&fit->ce_t, dim, 0);
	vector_init(&fit->be, 0);
	vector_init(&fit->duals, 0);
	vector_init(&fit->resid, dim);
	matrix_init(&fit->kkt, dim, dim);
	symeig_init(&fit->eig, dim, EIG_VEC);

	struct linesearch_ctrl lsctrl = LINESEARCH_CTRL0;
	fit->lsctrl = lsctrl;
	
	preprocess(fit);
}

void recv_fit_deinit(struct recv_fit *fit)
{
	assert(fit);
	
	symeig_deinit(&fit->eig);
	matrix_deinit(&fit->kkt);
	vector_deinit(&fit->resid);
	vector_deinit(&fit->duals);
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
