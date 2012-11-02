#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "xalloc.h"

#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"

#include "linesearch.h"
#include "recv_fit.h"

static size_t constr_count0(const struct constr *c);


size_t constr_count0(const struct constr *c)
{
	return c ? constr_count(c) : 0;
}

void recv_params_init(struct recv_params *x, const struct frame *f,
		      const struct constr *c)
{
	size_t coefs_dim = recv_coefs_dim(f);
	size_t nc = constr_count0(c);
	size_t dim = coefs_dim + nc;
	double *data = xmalloc(dim * sizeof(data[0]));
	recv_params_init_view(x, f, c, data);
	x->owner = 1;
}

 void recv_params_init_view(struct recv_params *x, const struct frame *f,
			    const struct constr *c, const double *data)
{
	x->all = (double *)data;
	x->dim = 0;

	recv_coefs_init_view(&x->coefs, f, x->all + x->dim);
	x->dim += x->coefs.dim;

	x->duals = x->all + x->dim;
	x->dim += constr_count0(c);

	x->owner = 0;
}


void recv_params_deinit(struct recv_params *x)
{
	if (x->owner)
		free(x->all);
}


static void resid_init(struct recv_fit_resid *resid, const struct frame *f,
		       const struct constr *c)
{
	recv_params_init(&resid->params, f, c);
	resid->norm2 = NAN;
}

static void resid_deinit(struct recv_fit_resid *resid)
{
	recv_params_deinit(&resid->params);
}


static void resid_set(struct recv_fit_resid *resid,
		      const struct recv_loglik *ll,
		      const struct constr *c,
		      const struct recv_params *params)
{
	const struct recv_model *m = recv_loglik_model(ll);
	size_t dim = recv_model_dim(m);
	size_t nc = constr_count0(c);
	const double *A = constr_all_wts(c);
	const double *b = constr_all_vals(c);

	memset(resid->params.all, 0, resid->params.dim * sizeof(double));
	if (!dim)
		return;

	/* r1 is the dual residual: grad(f) + A * nu,
	 *   where grad(f) = grad(nll)
	 *                 = -(score)
	 */
	double df = (double)recv_loglik_count(ll);
	recv_loglik_axpy_score(-1.0/df, ll, &resid->params.coefs);
	blas_dgemv(BLAS_NOTRANS, dim, nc, 1.0, A, dim,
		   params->duals, 1, 1.0, resid->params.coefs.all, 1);

	// r2 is the primal residual: A' * x - b
	blas_dcopy(nc, b, 1, resid->params.duals, 1);
	blas_dgemv(BLAS_TRANS, dim, nc, 1.0, A, dim,
		   params->coefs.all, 1, -1.0, resid->params.duals, 1);

	// compute ||r||^2
	double norm = blas_dnrm2(dim + nc, resid->params.all, 1);
	resid->norm2 = norm * norm;
}


static void eval_init(struct recv_fit_eval *eval, struct recv_model *m,
		      const struct constr *c)
{
	const struct frame *f = recv_model_frame(m);

	recv_loglik_init(&eval->loglik, m);
	recv_params_init(&eval->params, f, c);
	resid_init(&eval->resid, f, c);
}

static void eval_deinit(struct recv_fit_eval *eval)
{
	resid_deinit(&eval->resid);
	recv_params_deinit(&eval->params);
	recv_loglik_deinit(&eval->loglik);
}

static void _eval_update(struct recv_fit_eval *eval,
			 const struct messages *xmsgs,
			 const struct messages *ymsgs,
			 struct frame *frame, struct recv_model *model,
			 int moments)
{
	// clear frame; set model coefs
	struct history *h = frame_history(frame);
	history_clear(h);
	recv_model_set_coefs(model, &eval->params.coefs);
	recv_model_set_moments(model, moments);

	// set loglik
	recv_loglik_clear(&eval->loglik);
	struct messages_iter xit = messages_iter_make(xmsgs);
	struct messages_iter yit = messages_iter_make(ymsgs);
	double xt = -INFINITY;
	double yt = -INFINITY;
	const struct message *msg;
	size_t i, n;

	if (messages_iter_advance(&xit)) {
		xt = MESSAGES_TIME(xit);
	} else {
		xt = INFINITY;
	}

	while (messages_iter_advance(&yit)) {
		yt = MESSAGES_TIME(yit);

		while (xt < yt) {
			history_advance(h, xt);

			n = MESSAGES_COUNT(xit);
			for (i = 0; i < n; i++) {
				msg = MESSAGES_VAL(xit, i);
				history_add(h, msg);
			}

			if (messages_iter_advance(&xit)) {
				xt = MESSAGES_TIME(xit);
			} else {
				xt = INFINITY;
			}
		}

		history_advance(h, yt);
		n = MESSAGES_COUNT(yit);
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(yit, i);
			recv_loglik_add(&eval->loglik, frame, msg);
		}
	}
}

static void eval_set(struct recv_fit_eval *eval,
		     const struct recv_params *params,
		     const struct messages *xmsgs,
		     const struct messages *ymsgs,
		     struct frame *frame, struct recv_model *model, int moments)
{
	if (params) {
		memcpy(eval->params.all, params->all, eval->params.dim * sizeof(double));
	} else {
		memset(eval->params.all, 0, eval->params.dim * sizeof(double));
	}

	_eval_update(eval, xmsgs, ymsgs, frame, model, moments);
}

static void eval_step(struct recv_fit_eval *eval,
		      const struct recv_params *params0, double scale,
		      const struct recv_params *dir, const struct messages *xmsgs,
		      const struct messages *ymsgs,
		      struct frame *frame, struct recv_model *model, int moments)
{
	blas_dcopy(eval->params.dim, params0->all, 1, eval->params.all, 1);
	blas_daxpy(eval->params.dim, scale, dir->all, 1, eval->params.all, 1);
	_eval_update(eval, xmsgs, ymsgs, frame, model, moments);
}


static void kkt_init(struct recv_fit_kkt *kkt, const struct frame *f,
		     const struct constr *c)
{
	size_t dim = recv_coefs_dim(f);
	size_t nc = constr_count0(c);
	size_t n = dim + nc;
	size_t lwork = lapack_dsysv_lwork(n);

	kkt->factored = 0;
	kkt->matrix = xmalloc(n * n * sizeof(kkt->matrix[0]));;
	kkt->imat_buf = xmalloc(dim * (dim + 1) / 2 * sizeof(*kkt->imat_buf));
	kkt->ldl_ipiv = xmalloc(n * sizeof(kkt->ldl_ipiv[0]));;
	kkt->ldl_work = xmalloc(lwork * sizeof(kkt->ldl_work[0]));;
	kkt->ldl_lwork = lwork;
}

static void kkt_deinit(struct recv_fit_kkt *kkt)
{
	free(kkt->ldl_work);
	free(kkt->ldl_ipiv);
	free(kkt->imat_buf);
	free(kkt->matrix);
}

static void kkt_set(struct recv_fit_kkt *kkt, const struct recv_loglik *ll,
		    const struct constr *c)
{
	const struct recv_model *m = recv_loglik_model(ll);
	size_t dim = recv_model_dim(m);
	size_t nc = constr_count0(c);
	size_t n = dim + nc;

	kkt->factored = 0;

	memset(kkt->matrix, 0, n * n * sizeof(kkt->matrix[0]));

	/* k11 */
	memset(kkt->imat_buf, 0, dim * (dim + 1) / 2 * sizeof(kkt->imat_buf[0]));
	double df = (double)recv_loglik_count(ll);
	recv_loglik_axpy_imat(1.0/df, ll, kkt->imat_buf);
	enum blas_uplo f77uplo = RECV_LOGLIK_UPLO == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER;

	packed_dsctr(f77uplo, dim, kkt->imat_buf, kkt->matrix, n);

	if (!nc)
		return;

	if (RECV_LOGLIK_UPLO == BLAS_UPPER) { /* k12 */
		matrix_dtrans(dim, nc, c->wts, dim, kkt->matrix + dim, n);
	} else { /* k21 */
		lapack_dlacpy(LA_COPY_ALL, dim, nc, c->wts, dim,
			      kkt->matrix + dim * n, n);
	}
}

static void search_init(struct recv_fit_search *search, const struct frame *f,
			const struct constr *c)
{
	recv_params_init(&search->params, f, c);
}

static void search_deinit(struct recv_fit_search *search)
{
	recv_params_deinit(&search->params);
}

static void search_set(struct recv_fit_search *search,
		       const struct recv_fit_resid *resid,
		       struct recv_fit_kkt *kkt,
		       size_t n)
{
	assert(!kkt->factored);

	if (!n) {
		kkt->factored = 1;
		return;
	}

	/* determine the search direction */
	double *s = search->params.all;
	blas_dcopy(n, resid->params.all, 1, s, 1);
	blas_dscal(n, -1.0, s, 1);

	enum blas_uplo uplo = MLOGIT_COV_UPLO == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER;
	ptrdiff_t info = lapack_dsysv(uplo, n, 1, kkt->matrix, n,
				      kkt->ldl_ipiv, s, n, kkt->ldl_work,
				      kkt->ldl_lwork);
	kkt->factored = 1;

	assert(info == 0);
	(void)info;		// compile warning;
}

static void rgrad_init(struct recv_fit_rgrad *rgrad, const struct frame *f, const struct constr *c)
{
	recv_params_init(&rgrad->params, f, c);
}

static void rgrad_deinit(struct recv_fit_rgrad *rgrad)
{
	recv_params_deinit(&rgrad->params);
}

static void rgrad_set(struct recv_fit_rgrad *rgrad, size_t dim, size_t nc,
		      const struct recv_fit_kkt *kkt,
		      const struct recv_fit_resid *resid)
{
	assert(!kkt->factored);

	size_t n = dim + nc;
	enum blas_uplo f77uplo = RECV_LOGLIK_UPLO == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER;

	if (!n)
		return;

	blas_dsymv(f77uplo, n, 2.0, kkt->matrix, n, resid->params.all, 1, 0.0, rgrad->params.all, 1);
}

static enum recv_fit_task primal_dual_step(struct recv_fit *fit)
{
	size_t dim = recv_model_dim(&fit->model);
	size_t nc = constr_count0(fit->constr);
	size_t n = dim + nc;

	assert(dim > 0);

	// determine initial value and gradient
	double f0 = fit->cur->resid.norm2;
	//double g0 = -2 * f0;

	// stop early if residual is below tolerance
	if (f0 <= fit->ctrl.gtol * fit->ctrl.gtol)
		return RECV_FIT_CONV;

	// swap prev and cur evals
	struct recv_fit_eval *tmp = fit->prev;
	fit->prev = fit->cur;
	fit->cur = tmp;

	// determine the search direction
	search_set(&fit->search, &fit->prev->resid, &fit->kkt, n);
	size_t ismax = blas_idamax(dim, fit->search.params.coefs.all, 1) - 1;
	double smax = fabs(fit->search.params.coefs.all[ismax]);

	// set up the linesearch control parameters
	struct linesearch_ctrl ctrl = fit->ctrl.ls;
	ctrl.stpmax = MIN(ctrl.stpmax, 100 / MAX(1.0, smax));
	ctrl.stpmin = MIN(ctrl.stpmin, 1e-12 * ctrl.stpmax);
	double stp0 = MIN(1.0, 4.0 / MAX(1.0, smax));
	double stp = stp0;

	fprintf(stderr, "> smax = %.8f   stpmax = %.8f   stp0 = %.8f\n", smax,
		ctrl.stpmax, stp0);

	// perform a linesearch to reduce the residual norm     
	enum linesearch_task task;
	size_t it = 0;
	fit->step = stp0;
	//linesearch_start(&fit->ls, stp0, f0, g0, &ctrl);

	fprintf(stderr, ">                   f0: %.8f\n", f0);
	//fprintf(stderr, ">                   f0: %.8f  g0: %.8f\n", f0, g0);

	do {
		// compute the new trial step length
		it++;
		fit->step = stp;
		stp0 = stp;
		stp = stp0 * 0.5;

		// evaluate the function at the new point
		eval_step(fit->cur, &fit->prev->params, fit->step, &fit->search.params,
			  fit->xmsgs, fit->ymsgs, fit->frame, &fit->model, 1);

		resid_set(&fit->cur->resid, &fit->cur->loglik, fit->constr, &fit->cur->params);
		if (!isfinite(fit->cur->resid.norm2)) {
			goto domain_error_f;
		}

		double f = fit->cur->resid.norm2;
		//double g = vector_dot(&fit->search.vector, &fit->rgrad.vector);

		//if (!isfinite(g)) {
		//      goto domain_error_g;
		//}
		//fprintf(stderr, ">  stp: %.8f   f: %.8f   g: %.8f\n", fit->step, f, g);
		fprintf(stderr, ">  stp: %.8f   f: %.8f\n", fit->step, f);
		// stop (linesearch converged) or get a new trial step
		if (sqrt(f) <= (1 - 0.01 * fit->step) * sqrt(f0)) {
			task = LINESEARCH_CONV;
		} else if (stp < fit->ctrl.xtol) {
			task = LINESEARCH_WARN_XTOL;
		} else {
			task = LINESEARCH_STEP;
		}
		//task = linesearch_advance(&fit->ls, f, g);
		continue;

domain_error_f:
		//domain_error_g:
		assert(0 && "DOMAIN ERROR");
		task = LINESEARCH_STEP;
	} while (it < fit->ctrl.ls_maxit && task == LINESEARCH_STEP);	// && !linesearch_sdec(&fit->ls));

	// compute the kkt matrix and residual gradient at the new point
	eval_set(fit->cur, &fit->cur->params, fit->xmsgs, fit->ymsgs, fit->frame, &fit->model, 2);
	kkt_set(&fit->kkt, &fit->cur->loglik, fit->constr);
	rgrad_set(&fit->rgrad, dim, nc, &fit->kkt, &fit->cur->resid); // TODO: remove this?

	if (task == LINESEARCH_WARN_XTOL) {
		return RECV_FIT_ERR_XTOL;
	} else if (task == LINESEARCH_CONV) {
		return RECV_FIT_STEP;
	} else {
		return RECV_FIT_ERR_LNSRCH;
	}
}

void recv_fit_init(struct recv_fit *fit,
		   struct frame *f,
		   struct constr *c,
		   const struct messages *xmsgs,
		   const struct messages *ymsgs,
		   const struct recv_fit_ctrl *ctrl)
{
	assert(fit);
	assert(xmsgs);
	assert(f);
	assert(!c || recv_coefs_dim(f) == constr_dim(c));
	assert(messages_max_to(xmsgs) < frame_send_count(f));
	assert(messages_max_from(xmsgs) < frame_recv_count(f));
	assert(!ymsgs || (messages_max_to(ymsgs) < frame_send_count(f)));
	assert(!ymsgs || (messages_max_from(ymsgs) < frame_recv_count(f)));
	assert(!ctrl || recv_fit_ctrl_valid(ctrl));

	struct history *h = frame_history(f);

	if (!ctrl) {
		fit->ctrl = RECV_FIT_CTRL0;
	} else {
		fit->ctrl = *ctrl;
	}

	history_clear(h);
	fit->frame = f;
	fit->xmsgs = xmsgs;
	fit->ymsgs = ymsgs ? ymsgs : xmsgs;

	recv_model_init(&fit->model, f, NULL);
	fit->constr = c;
	eval_init(&fit->eval[0], &fit->model, c);
	eval_init(&fit->eval[1], &fit->model, c);
	fit->prev = &fit->eval[0];
	fit->cur = &fit->eval[1];
	kkt_init(&fit->kkt, f, c);
	search_init(&fit->search, f, c);
	rgrad_init(&fit->rgrad, f, c);

	// evaluate the model at zero
	eval_set(fit->cur, NULL, fit->xmsgs, fit->ymsgs, fit->frame, &fit->model, 2);

	fit->dev0 = recv_fit_dev(fit);
}


enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct recv_params *params0)
{
	if (params0 != NULL) {
		eval_set(fit->cur, params0, fit->xmsgs, fit->ymsgs, fit->frame,
			 &fit->model, 2);
	}

	resid_set(&fit->cur->resid, &fit->cur->loglik, fit->constr, &fit->cur->params);
	kkt_set(&fit->kkt, &fit->cur->loglik, fit->constr);
	fit->step = NAN;

	// determine initial value and gradient
	double f0 = fit->cur->resid.norm2;
	//double g0 = -2 * f0;

	// int in_domain = isfinite(f0);


	// stop early if residual is below tolerance
	if (f0 <= fit->ctrl.gtol * fit->ctrl.gtol) {
		fit->task = RECV_FIT_CONV;
	} else {
		fit->task = RECV_FIT_STEP;
	}

	return fit->task;
}

void recv_fit_deinit(struct recv_fit *fit)
{
	assert(fit);

	rgrad_deinit(&fit->rgrad);
	search_deinit(&fit->search);
	kkt_deinit(&fit->kkt);
	eval_deinit(&fit->eval[1]);
	eval_deinit(&fit->eval[0]);
	recv_model_deinit(&fit->model);
}


/*
size_t recv_fit_add_constr_identify(struct recv_fit *fit)
{
	const struct recv_loglik *ll = &fit->cur->loglik;
	const struct recv_model *m = recv_loglik_model(ll);
	size_t n = recv_model_dim(m);
	double *imatp = xmalloc(n * (n + 1) / 2 * sizeof(imatp[0]));
	memset(imatp, 0, n * (n + 1) / 2 * sizeof(imatp[0]));
	recv_loglik_axpy_imat(1.0, ll, imatp);

	size_t nadd = constr_add_identify(fit->constr, imatp, RECV_LOGLIK_UPLO);
	free(imatp);
	return nadd;
}
*/


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

const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit)
{
	assert(fit);
	return &fit->cur->loglik;
}

double recv_fit_dev(const struct recv_fit *fit)
{
	assert(fit);

	const struct recv_loglik *ll = &fit->cur->loglik;
	double dev = recv_loglik_dev(ll);
	return dev;
}

const struct recv_coefs *recv_fit_coefs(const struct recv_fit *fit)
{
	assert(fit);
	return &fit->cur->params.coefs;
}

const double *recv_fit_duals(const struct recv_fit *fit)
{
	assert(fit);
	return fit->cur->params.duals;
}

double recv_fit_step(const struct recv_fit *fit)
{
	assert(fit);
	return fit->step;
}

double recv_fit_grad_norm2(const struct recv_fit *fit)
{
	assert(fit);
	return fit->cur->resid.norm2;
}

double recv_fit_dev0(const struct recv_fit *fit)
{
	return fit->dev0;
}
