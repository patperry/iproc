#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linesearch.h"
#include "util.h"
#include "recv_fit.h"

// sets model and loglik
// advances frame to end of messages
// assumes coefs has already been set
static void evaluate_loglik(struct recv_fit *fit)
{
	frame_clear(&fit->frame);

	// update coefs
	recv_model_set_coefs(&fit->model, &fit->coefs);

	// update loglik
	recv_loglik_clear(&fit->loglik);
	recv_loglik_add_all(&fit->loglik, &fit->frame, fit->msgs);
}

// sets scale
// assumes evaluate has already been called
static void cohort_compute_scale(struct recv_fit *fit, ssize_t c)
{
	struct recv_fit_cohort *cohort = &fit->cohorts[c];	
	const double vartol = fit->ctrl.vartol;
	ssize_t i, dim = matrix_nrow(&fit->coefs);

	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik, c);
	const struct matrix *imat = &info->imat;

	/* get the empirical variances of the covariates */
	assert(vector_dim(&cohort->scale) == matrix_diag_dim(imat, 0));
	matrix_get_diag(imat, 0, vector_to_ptr(&cohort->scale));

	/* define the variance to be 1 when the variance is 0 */
	for (i = 0; i < dim; i++) {
		double *s = vector_item_ptr(&cohort->scale, i);
		if (*s < vartol)
			*s = 1;
	}

	/* otherwise, take the square root of the variance */
	vector_sqrt(&cohort->scale);
}

static void compute_scale(struct recv_fit *fit)
{
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	for (c = 0; c < n; c++) {
		cohort_compute_scale(fit, c);
	}
}

// sets ce_t, be, and ne
// reinits kkt and resid, duals, and ldlfac
// assumes scale is initalized and evaluate has already been called
static void cohort_compute_constraints(struct recv_fit *fit, ssize_t c)
{
	struct recv_fit_cohort *cohort = &fit->cohorts[c];
	const double eigtol = fit->ctrl.eigtol;
	ssize_t i, dim = matrix_nrow(&fit->coefs);
	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik, c);
	const struct matrix *imat = &info->imat;

	/* compute the eigendecomp of the empirical covariance matrix;
	 * kkt matrix and resid vector are un-initialized/free, so we
	 * can use them for temporary storage.
	 */
	struct matrix evec = matrix_slice(&cohort->kkt, 0, 0, dim, dim);
	struct vector eval = vector_slice(&cohort->resid, 0, dim);

	matrix_assign_copy(&evec, TRANS_NOTRANS, imat);

	/* compute its eigendecomposition */
	assert(symeig_dim(&fit->eig) == dim);
	assert(symeig_job(&fit->eig) == EIG_VEC);
	bool ok = symeig_factor(&fit->eig, UPLO_LOWER, &evec, &eval);
	assert(ok);

	/* compute its rank */
	for (i = 0; i < dim; i++) {
		double l = vector_item(&eval, i);
		if (l > eigtol)
			break;
	}
	// rank = dim - i;

	/* define the constraints */
	struct matrix u = matrix_slice_cols(&evec, 0, i);
	matrix_reinit(&cohort->ce, dim, i);
	matrix_assign_copy(&cohort->ce, TRANS_NOTRANS, &u);
	vector_reinit(&cohort->be, i);
	vector_fill(&cohort->be, 0);
	cohort->ne = i;

	/* update dims of kkt matrix, residuals and dual parameters */
	matrix_reinit(&cohort->kkt, dim + i, dim + i);
	vector_reinit(&cohort->resid, dim + i);
	vector_reinit(&cohort->params0, dim + i);
	vector_reinit(&cohort->params, dim + i);	
	vector_reinit(&cohort->search, dim + i);
	vector_reinit(&cohort->grad_rss, dim + i);	

	/* update ldlfac dimension */
	ldlfac_reinit(&cohort->ldl, dim + i);
}


static void compute_constraints(struct recv_fit *fit)
{
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	for (c = 0; c < n; c++) {
		cohort_compute_constraints(fit, c);
	}
}


// sets resid
// assumes evaluate and preprocess have already been called
static void cohort_compute_resid(struct recv_fit *fit, ssize_t c)
{
	struct recv_fit_cohort *cohort = &fit->cohorts[c];
	ssize_t dim = recv_model_dim(&fit->model);
	ssize_t ne = cohort->ne;

	assert(vector_dim(&cohort->be) == ne);

	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik, c);
	const struct vector *score = &info->score;

	/* r1 is the dual residual: grad(f) + ce * nu,
	 * where grad(f) = grad(nll) + lambda * x
	 *               = -[score - lambda * x]
	 */
	struct vector r1 = vector_slice(&cohort->resid, 0, dim);
	struct vector primals = vector_slice(&cohort->params, 0, dim);
	struct vector duals = vector_slice(&cohort->params, dim, ne);
	vector_assign_copy(&r1, score);
	vector_axpy(-fit->penalty, &primals, &r1);
	matrix_mul(1.0, TRANS_NOTRANS, &cohort->ce, &duals, -1.0, &r1);

	/* r2 is the primal residual: ce * x - be */
	struct vector r2 = vector_slice(&cohort->resid, dim, ne);
	vector_assign_copy(&r2, &cohort->be);
	matrix_mul(1.0, TRANS_TRANS, &cohort->ce, &primals, -1.0, &r2);
}

static void compute_resid(struct recv_fit *fit)
{
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	for (c = 0; c < n; c++) {
		cohort_compute_resid(fit, c);
	}
}

// sets kkt
// assumes evaluate and preproeces have already been called
static void cohort_compute_kkt(struct recv_fit *fit, ssize_t c)
{
	struct recv_fit_cohort *cohort = &fit->cohorts[c];
	ssize_t dim = recv_model_dim(&fit->model);
	ssize_t ne = cohort->ne;

	assert(matrix_nrow(&cohort->kkt) == dim + ne);
	assert(matrix_ncol(&cohort->kkt) == dim + ne);

	const struct recv_loglik_info *info = recv_loglik_info(&fit->loglik, c);
	const struct matrix *imat = &info->imat;

	/* k11 is the hessian */
	struct matrix k11 = matrix_slice(&cohort->kkt, 0, 0, dim, dim);
	matrix_assign_copy(&k11, TRANS_NOTRANS, imat);
	ssize_t i;
	for (i = 0; i < dim; i++) {
		*matrix_item_ptr(&k11, i, i) += fit->penalty;
	}

	/* k12 is the transpose of the equality constraint matrix */
	struct matrix k12 = matrix_slice(&cohort->kkt, 0, dim, dim, ne);
	matrix_assign_copy(&k12, TRANS_NOTRANS, &cohort->ce);

	/* k22 is zero */
	struct matrix k22 = matrix_slice(&cohort->kkt, dim, dim, ne, ne);
	matrix_fill(&k22, 0.0);

	/* k21 is never referenced */
#ifndef NDEBUG
	struct matrix k21 = matrix_slice(&cohort->kkt, dim, 0, ne, dim);
	matrix_fill(&k21, 0);
#endif
}

static void compute_kkt(struct recv_fit *fit)
{
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	for (c = 0; c < n; c++) {
		cohort_compute_kkt(fit, c);
	}
}

static double resid_norm2(const struct recv_fit *fit)
{
	double rss = 0.0;
	
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	for (c = 0; c < n; c++) {
		rss += vector_norm2(&fit->cohorts[c].resid);
	}
	return rss;
}


// sets scale, ce, be
// reinits kkt and nresid
// destroys model and coefs 
static void preprocess(struct recv_fit *fit, const struct matrix *coefs0)
{
	ssize_t c, nc = recv_model_cohort_count(&fit->model);
	/* compute the empirical covariance matrix of the covariates
	 * by using the null model (coefs = 0, all receivers equally
	 * likely
	 */
	matrix_fill(&fit->coefs, 0);
	evaluate_loglik(fit);

	const struct recv_loglik *ll = recv_fit_loglik(fit);
	
	for (c = 0; c < nc; c++) {
		const struct recv_loglik_info *info = recv_loglik_info(ll, c);
		matrix_assign_copy(&fit->cohorts[c].imat0, TRANS_NOTRANS, &info->imat);
		vector_assign_copy(&fit->cohorts[c].score0, &info->score);
		fit->cohorts[c].dev0 = info->dev;
		
		cohort_compute_scale(fit, c);
		cohort_compute_constraints(fit, c);

	}

	if (coefs0) {
		matrix_assign_copy(&fit->coefs, TRANS_NOTRANS, coefs0);
		evaluate_loglik(fit);
	} else {
		assert(vector_norm(&fit->coefs.data) == 0);
	}
	
	for (c = 0; c < nc; c++) {
		ssize_t dim = recv_model_dim(&fit->model);
		ssize_t ne = fit->cohorts[c].ne;
		struct vector primals = vector_slice(&fit->cohorts[c].params, 0, dim);
		struct vector duals = vector_slice(&fit->cohorts[c].params, dim, ne);

		struct vector coefs_c = matrix_col(&fit->coefs, c);
		vector_assign_copy(&primals, &coefs_c);
		vector_fill(&duals, 0.0);
	}

	compute_resid(fit);
	compute_kkt(fit);
	fit->task = RECV_FIT_STEP;
	fit->rss = resid_norm2(fit);
}

static void cohort_compute_search(struct recv_fit *fit, ssize_t c)
{
	ssize_t dim = recv_model_dim(&fit->model);
	ssize_t ne = fit->cohorts[c].ne;

	/* determine the search direction */
	vector_assign_copy(&fit->cohorts[c].search, &fit->cohorts[c].resid);
	struct matrix smat = matrix_make(&fit->cohorts[c].search, dim + ne, 1);
	ssize_t info = ldlfac_solve(&fit->cohorts[c].ldl, UPLO_UPPER,
				    &fit->cohorts[c].kkt, &smat);
	assert(info == 0);
	vector_scale(&fit->cohorts[c].search, -1.0);
	
	/* update the initial position */
	vector_assign_copy(&fit->cohorts[c].params0, &fit->cohorts[c].params);
	
}

static void compute_search(struct recv_fit *fit)
{
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	for (c = 0; c < n; c++) {
		cohort_compute_search(fit, c);
	}
}


static void take_step(struct recv_fit *fit, double stp, double *rssp, double *grssp)
{
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	ssize_t dim = recv_model_dim(&fit->model);
	
	/* take a step to get new primal and dual variables (and coefficients) */
	for (c = 0; c < n; c++) {
		vector_assign_copy(&fit->cohorts[c].params, &fit->cohorts[c].params0);
		vector_axpy(stp, &fit->cohorts[c].search, &fit->cohorts[c].params);

		struct vector primals = vector_slice(&fit->cohorts[c].params, 0, dim);
		struct vector coefs_c = matrix_col(&fit->coefs, c);
		vector_assign_copy(&coefs_c, &primals);
	}
	
	/* compute the fit, residuals, and kkt matrix at the new coefficients */
	evaluate_loglik(fit);
	compute_resid(fit);
	compute_kkt(fit);

	/* compute the gradient of the squared residual norm */
	double rss = 0.0, grss = 0.0;
	
	for (c = 0; c < n; c++) {
		ssize_t ne = fit->cohorts[c].ne;
		struct vector cgrad = vector_slice(&fit->cohorts[c].grad_rss, 0, dim);
		struct vector dgrad = vector_slice(&fit->cohorts[c].grad_rss, dim, ne);	
		struct vector cresid = vector_slice(&fit->cohorts[c].resid, 0, dim);
		struct vector dresid = vector_slice(&fit->cohorts[c].resid, dim, ne);
		struct matrix k11 = matrix_slice(&fit->cohorts[c].kkt, 0, 0, dim, dim);
	
		matrix_mul(1.0, TRANS_NOTRANS, &k11, &cresid, 0.0, &cgrad);
		matrix_mul(1.0, TRANS_NOTRANS, &fit->cohorts[c].ce, &dresid, 1.0, &cgrad);
		matrix_mul(1.0, TRANS_TRANS, &fit->cohorts[c].ce, &cresid, 0.0, &dgrad);
		vector_scale(&cgrad, 2.0);
		vector_scale(&dgrad, 2.0);
		
		/* evaluate the new squared residual norm and directional derivative */
		fit->cohorts[c].rss = vector_norm2(&fit->cohorts[c].resid);
		fit->cohorts[c].grss = vector_dot(&fit->cohorts[c].search, &fit->cohorts[c].grad_rss);
		
		rss += fit->cohorts[c].rss;
		grss += fit->cohorts[c].grss;
	}
	
	assert(isfinite(rss));
	assert(isfinite(grss));
	
	*rssp = rss;
	*grssp = grss;
}

static enum recv_fit_task primal_dual_step(struct recv_fit *fit)
{

	double f0 = fit->rss;
	double g0 = -2 * f0;

	assert(isfinite(f0));
	assert(isfinite(g0));

	if (-g0 < fit->ctrl.gtol * fit->ctrl.gtol)
		return RECV_FIT_CONV;	/* residual is 0 */

	compute_search(fit);

	/* perform a linesearch to reduce the residual norm */
	enum linesearch_task task;
	ssize_t it = 0;

	//printf("f0: %.22f g0: %.22f\n", f0, g0);


	
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	double xmax = 0.0;
	for (c = 0; c < n; c++) {
		xmax = MAX(xmax, vector_max_abs(&fit->cohorts[c].search));
	}
	
	struct linesearch_ctrl ctrl = fit->ctrl.ls;	
	ctrl.stpmax = MIN(ctrl.stpmax, 1000 / MAX(1.0, xmax));
	ctrl.stpmin = MIN(ctrl.stpmin, 1e-12 * ctrl.stpmax);
	double stp0 = MIN(1.0, ctrl.stpmin + 0.5 * (ctrl.stpmax - ctrl.stpmin));
	linesearch_start(&fit->ls, stp0, f0, g0, &ctrl);

	do {
		it++;
		fit->step = linesearch_step(&fit->ls);
		take_step(fit, fit->step, &fit->rss, &fit->grss);

		task = linesearch_advance(&fit->ls, fit->rss, fit->grss);

		//printf("f: %.22f g: %.22f  stp: %.22f\n", f, g, stp);
	} while (it < fit->ctrl.ls_maxit
		 && task == LINESEARCH_STEP && !linesearch_sdec(&fit->ls));

	if (fit->step < fit->ctrl.xtol) {
		return RECV_FIT_ERR_XTOL;
	} else if (linesearch_sdec(&fit->ls)) {
		return RECV_FIT_STEP;
	} else {
		return RECV_FIT_ERR_LNSRCH;
	}
}

static void cohort_init(struct recv_fit_cohort *cohort, ssize_t dim)
{
	matrix_init(&cohort->imat0, dim, dim);
	vector_init(&cohort->score0, dim);
	vector_init(&cohort->scale, dim);
	matrix_init(&cohort->ce, dim, 0);
	vector_init(&cohort->be, 0);
	vector_init(&cohort->params0, dim);
	vector_init(&cohort->params, dim);	
	vector_init(&cohort->resid, dim);
	matrix_init(&cohort->kkt, dim, dim);
	vector_init(&cohort->search, dim);
	vector_init(&cohort->grad_rss, dim);	
	ldlfac_init(&cohort->ldl, dim);
}

void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct actors *senders,		   
		   const struct matrix *coefs0,
		   const struct recv_fit_ctrl *ctrl)
{
	assert(fit);
	assert(msgs);
	assert(design);
	assert(messages_max_to(msgs) < design_send_count(design));
	assert(messages_max_from(msgs) < design_recv_count(design));
	assert(!coefs0 || matrix_nrow(coefs0) == design_recv_dim(design));
	assert(!ctrl || recv_fit_ctrl_valid(ctrl));

	ssize_t dim = design_recv_dim(design);
	ssize_t ncohort = actors_cohort_count(senders);

	if (!ctrl) {
		fit->ctrl = RECV_FIT_CTRL0;
	} else {
		fit->ctrl = *ctrl;
	}

	fit->design = design;
	fit->senders = senders;
	fit->msgs = msgs;
	fit->penalty = 0.0;

	frame_init(&fit->frame, design);
	matrix_init(&fit->coefs, dim, ncohort);
	recv_model_init(&fit->model, &fit->frame, fit->senders, &fit->coefs);
	recv_loglik_init(&fit->loglik, &fit->model);
	fit->step = NAN;
	symeig_init(&fit->eig, dim, EIG_VEC);

	ssize_t c;
	struct recv_fit_cohort *cohorts = xcalloc(ncohort, sizeof(*cohorts));
	for (c = 0; c < ncohort; c++) {
		cohort_init(&cohorts[c], dim);
	}
	fit->cohorts = cohorts;
	
	preprocess(fit, coefs0);
}

static void cohort_deinit(struct recv_fit_cohort *cohort)
{
	ldlfac_deinit(&cohort->ldl);
	vector_deinit(&cohort->grad_rss);	
	vector_deinit(&cohort->search);
	matrix_deinit(&cohort->kkt);
	vector_deinit(&cohort->resid);
	vector_deinit(&cohort->params);
	vector_deinit(&cohort->params0);	
	vector_deinit(&cohort->be);
	matrix_deinit(&cohort->ce);
	vector_deinit(&cohort->scale);
	vector_deinit(&cohort->score0);
	matrix_deinit(&cohort->imat0);
}

void recv_fit_deinit(struct recv_fit *fit)
{
	assert(fit);

	
	ssize_t c, n = recv_model_cohort_count(&fit->model);
	struct recv_fit_cohort *cohorts = fit->cohorts;
	for (c = 0; c < n; c++) {
		cohort_deinit(&cohorts[c]);
	}
	xfree(cohorts);

	symeig_deinit(&fit->eig);
	recv_loglik_deinit(&fit->loglik);
	recv_model_deinit(&fit->model);
	matrix_deinit(&fit->coefs);
	frame_deinit(&fit->frame);
}

enum recv_fit_task recv_fit_advance(struct recv_fit *fit)
{
	assert(fit);
	assert(fit->task == RECV_FIT_STEP);

	fit->task = primal_dual_step(fit);
	return fit->task;
}

const char *recv_fit_errmsg(const struct recv_fit *fit)
{
	assert(fit);

	switch (fit->task) {
	case RECV_FIT_STEP:
		return "optimization in progress";
	case RECV_FIT_ERR_LNSRCH:
		return "linesearch failed";
	case RECV_FIT_ERR_XTOL:
		return "step size is less than tolerance";
	case RECV_FIT_CONV:
		return NULL;
	}

	assert(0);
	return NULL;
}

ssize_t recv_fit_rank(const struct recv_fit *fit, ssize_t c)
{
	assert(fit);
	ssize_t dim = matrix_nrow(&fit->coefs);
	ssize_t ne = vector_dim(&fit->cohorts[c].be);

	return dim - ne;
}

const struct matrix *recv_fit_ce(const struct recv_fit *fit, ssize_t c)
{
	assert(fit);
	return &fit->cohorts[c].ce;
}

const struct vector *recv_fit_be(const struct recv_fit *fit, ssize_t c)
{
	assert(fit);
	return &fit->cohorts[c].be;
}

double recv_fit_dev0(const struct recv_fit *fit, ssize_t c)
{
	assert(fit);
	return fit->cohorts[c].dev0;
}

const struct vector *recv_fit_score0(const struct recv_fit *fit, ssize_t c)
{
	assert(fit);
	return &fit->cohorts[c].score0;
}

const struct matrix *recv_fit_imat0(const struct recv_fit *fit, ssize_t c)
{
	assert(fit);
	return &fit->cohorts[c].imat0;
}

const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit)
{
	assert(fit);
	return &fit->loglik;
}

const struct matrix *recv_fit_coefs(const struct recv_fit *fit)
{
	assert(fit);
	return &fit->coefs;
}

double recv_fit_step(const struct recv_fit *fit)
{
	assert(fit);
	return fit->step;
}

double recv_fit_grad_norm2(const struct recv_fit *fit)
{
	assert(fit);
	return fit->rss;
}
