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

static void params_init(struct recv_fit_params *x, const struct frame *f, size_t ne);
static void params_init_view(struct recv_fit_params *x, const struct frame *f, size_t ne, const double *data);
static void params_reinit(struct recv_fit_params *x, const struct frame *f, size_t ne);
static void params_clear(struct recv_fit_params *x);
static void params_deinit(struct recv_fit_params *x);

static void params_init(struct recv_fit_params *x, const struct frame *f, size_t ne)
{
	size_t coefs_dim = recv_coefs_dim(f);
	size_t dim = coefs_dim + ne;
	double *data = xmalloc(dim * sizeof(data[0]));
	params_init_view(x, f, ne, data);
	x->owner = 1;
	params_clear(x);
}

static void params_init_view(struct recv_fit_params *x, const struct frame *f, size_t ne, const double *data)
{
	x->all = (double *)data;
	x->dim = 0;

	recv_coefs_init_view(&x->coefs, f, x->all + x->dim);
	x->dim += x->coefs.dim;

	x->duals = x->all + x->dim;
	x->dim += ne;

	x->owner = 0;
}


static void params_reinit(struct recv_fit_params *x, const struct frame *f, size_t ne)
{
	assert(x->owner);

	size_t dim0 = x->coefs.dim;
	size_t ne0 = x->dim - dim0;
	size_t dim = dim0 + ne;

	double *data = xrealloc(x->all, dim * sizeof(data[0]));
	params_init_view(x, f, ne, data);
	x->owner = 1;

	if (ne > ne0) {
		memset(x->duals + ne0, 0, (ne - ne0) * sizeof(*x->duals));
	}
}


static void params_clear(struct recv_fit_params *x)
{
	memset(x->all, 0, x->dim * sizeof(*x->all));
}


static void params_deinit(struct recv_fit_params *x)
{
	if (x->owner)
		free(x->all);
}


void recv_fit_params_init(struct recv_fit_params *params, const struct recv_fit *fit)
{
	const struct frame *f = fit->frame;
	size_t ne = recv_fit_constr_count(fit);
	params_init(params, f, ne);
}


void recv_fit_params_deinit(struct recv_fit_params *params)
{
	params_deinit(params);
}


static void constr_init(struct recv_fit_constr *c)
{
	c->wts = NULL;
	c->vals = NULL;
	c->n = 0;
	c->nmax = 0;
}


static void constr_deinit(struct recv_fit_constr *c)
{
	free(c->vals);
	free(c->wts);
}

static void constr_grow(struct recv_fit_constr *c, size_t dim)
{
	if (c->n == c->nmax) {
		size_t nmax = ARRAY_GROW1(c->nmax, SIZE_MAX);
		c->wts = xrealloc(c->wts, nmax * dim * sizeof(*c->wts));
		c->vals = xrealloc(c->vals, nmax * sizeof(*c->vals));
		c->nmax = nmax;
	}
}

static size_t constr_count(const struct recv_fit_constr *c)
{
	return c->n;
}

static void constr_add(struct recv_fit_constr *c, size_t dim, const double *wts,
		       double val)
{
	constr_grow(c, dim);

	size_t n = c->n;
	memcpy(c->wts + n * dim, wts, dim * sizeof(*c->wts));
	c->vals[n] = val;
	c->n = n + 1;
}

static void constr_add_set(struct recv_fit_constr *c, size_t dim, size_t i, double val)
{
	assert(c);
	assert(i < dim);
	assert(isfinite(val));

	constr_grow(c, dim);

	size_t n = c->n;
	memset(c->wts + n * dim, 0, dim * sizeof(*c->wts));
	c->wts[i + n * dim] = 1.0;
	c->vals[n] = val;
	c->n = n + 1;
}

static void constr_add_eq(struct recv_fit_constr *c, size_t dim, size_t i1, size_t i2)
{
	assert(c);
	assert(i1 != i2);
	assert(i1 < dim);
	assert(i2 < dim);

	constr_grow(c, dim);

	size_t n = c->n;
	memset(c->wts + n * dim, 0, dim * sizeof(*c->wts));
	c->wts[i1 + n * dim] = +1.0;
	c->wts[i2 + n * dim] = -1.0;
	c->vals[n] = 0.0;
	c->n = n + 1;
}

static void constr_get(const struct recv_fit_constr *c, size_t dim, size_t i,
		       const double **weightsp, double *valp)
{
	assert(i < constr_count(c));

	const double *wts = c->wts + i * dim;
	double val = c->vals[i];

	if (weightsp)
		*weightsp = wts;
	if (valp)
		*valp = val;
}

static void resid_init(struct recv_fit_resid *resid, const struct frame *f,
		       size_t ne)
{
	params_init(&resid->params, f, ne);
	resid->norm2 = NAN;
}

static void resid_deinit(struct recv_fit_resid *resid)
{
	params_deinit(&resid->params);
}

static void resid_reinit(struct recv_fit_resid *resid, const struct frame *f, size_t ne)
{
	params_reinit(&resid->params, f, ne);
	resid->norm2 = NAN;
}


static void resid_set(struct recv_fit_resid *resid,
		      const struct recv_loglik *ll,
		      const struct recv_fit_constr *ce,
		      const struct recv_fit_params *params)
{
	const struct recv_model *m = recv_loglik_model(ll);
	size_t dim = recv_model_dim(m);
	size_t ne = constr_count(ce);

	params_clear(&resid->params);
	if (!dim)
		return;

	/* r1 is the dual residual: grad(f) + ce * nu,
	 *   where grad(f) = grad(nll)
	 *                 = -(score)
	 */
	double df = (double)recv_loglik_count(ll);
	recv_loglik_axpy_score(-1.0/df, ll, &resid->params.coefs);
	blas_dgemv(BLAS_NOTRANS, dim, ne, 1.0, ce->wts, dim,
		   params->duals, 1, 1.0, resid->params.coefs.all, 1);

	// r2 is the primal residual: ce' * x - be
	blas_dcopy(ne, ce->vals, 1, resid->params.duals, 1);
	blas_dgemv(BLAS_TRANS, dim, ne, 1.0, ce->wts, dim,
		   params->coefs.all, 1, -1.0, resid->params.duals, 1);

	// compute ||r||^2
	double norm = blas_dnrm2(dim + ne, resid->params.all, 1);
	resid->norm2 = norm * norm;
}


static void eval_init(struct recv_fit_eval *eval, struct recv_model *m,
		      size_t ne)
{
	const struct frame *f = recv_model_frame(m);

	recv_loglik_init(&eval->loglik, m);
	params_init(&eval->params, f, ne);
	resid_init(&eval->resid, f, ne);
}

static void eval_deinit(struct recv_fit_eval *eval)
{
	resid_deinit(&eval->resid);
	params_deinit(&eval->params);
	recv_loglik_deinit(&eval->loglik);
}

static void eval_reinit(struct recv_fit_eval *eval, const struct frame *f,
		        size_t ne)
{
	params_reinit(&eval->params, f, ne);
	resid_reinit(&eval->resid, f, ne);
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
		     const struct recv_fit_params *params,
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
		      const struct recv_fit_params *params0, double scale,
		      const struct recv_fit_params *dir, const struct messages *xmsgs,
		      const struct messages *ymsgs,
		      struct frame *frame, struct recv_model *model, int moments)
{
	blas_dcopy(eval->params.dim, params0->all, 1, eval->params.all, 1);
	blas_daxpy(eval->params.dim, scale, dir->all, 1, eval->params.all, 1);
	_eval_update(eval, xmsgs, ymsgs, frame, model, moments);
}

static void kkt_reinit(struct recv_fit_kkt *kkt, const struct frame *f, size_t ne)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dim = design_dim(r) + design2_dim(d);
	size_t n = dim + ne;
	size_t lwork = lapack_dsysv_lwork(n);

	kkt->factored = 0;
	kkt->matrix = xrealloc(kkt->matrix, n * n * sizeof(kkt->matrix[0]));
	kkt->ldl_ipiv = xrealloc(kkt->ldl_ipiv, n * sizeof(kkt->ldl_ipiv[0]));
	kkt->ldl_work = xrealloc(kkt->ldl_work, lwork * sizeof(kkt->ldl_work[0]));
	kkt->ldl_lwork = lwork;
}

static void kkt_init(struct recv_fit_kkt *kkt, const struct frame *f, size_t ne)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dim = design_dim(r) + design2_dim(d);

	kkt->factored = 0;
	kkt->matrix = NULL;
	kkt->imat_buf = xmalloc(dim * (dim + 1) / 2 * sizeof(*kkt->imat_buf));
	kkt->ldl_ipiv = NULL;
	kkt->ldl_work = NULL;
	kkt->ldl_lwork = 0;
	kkt_reinit(kkt, f, ne);
}

static void kkt_deinit(struct recv_fit_kkt *kkt)
{
	free(kkt->ldl_work);
	free(kkt->ldl_ipiv);
	free(kkt->imat_buf);
	free(kkt->matrix);
}

static void kkt_set(struct recv_fit_kkt *kkt, const struct recv_loglik *ll,
		    const struct recv_fit_constr *ce)
{
	const struct recv_model *m = recv_loglik_model(ll);
	size_t i, dim = recv_model_dim(m);
	size_t ne = constr_count(ce);
	size_t n = dim + ne;

	memset(kkt->matrix, 0, n * n * sizeof(kkt->matrix[0]));

	// k11
	memset(kkt->imat_buf, 0, dim * (dim + 1) / 2 * sizeof(kkt->imat_buf[0]));
	double df = (double)recv_loglik_count(ll);
	recv_loglik_axpy_imat(1.0/df, ll, kkt->imat_buf);
	double *dst = kkt->matrix;
	const double *src = kkt->imat_buf;

	if (MLOGIT_COV_UPLO == BLAS_UPPER) {
		for (i = 0; i < dim; i++) {
			size_t rowlen = dim - i;
			
			memcpy(dst + i, src, rowlen * sizeof(dst[0]));
			dst += n;
			src += rowlen;
		}

		// k12
		if (ne)
			matrix_dtrans(dim, ne, ce->wts, dim, kkt->matrix + dim, n);
	} else {
		for (i = 0; i < dim; i++) {
			size_t rowlen = i + 1;

			memcpy(dst, src, rowlen * sizeof(dst[0]));
			dst += n;
			src += rowlen;
		}

		// k21
		if (ne)
			lapack_dlacpy(LA_COPY_ALL, dim, ne, ce->wts, dim, kkt->matrix + dim * n, n);

	}

	kkt->factored = 0;
}

static void search_init(struct recv_fit_search *search, const struct frame *f, size_t ne)
{
	params_init(&search->params, f, ne);
}

static void search_reinit(struct recv_fit_search *search, const struct frame *f, size_t ne)
{
	params_reinit(&search->params, f, ne);
}

static void search_deinit(struct recv_fit_search *search)
{
	params_deinit(&search->params);
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

static void rgrad_init(struct recv_fit_rgrad *rgrad, const struct frame *f, size_t ne)
{
	params_init(&rgrad->params, f, ne);
}

static void rgrad_reinit(struct recv_fit_rgrad *rgrad, const struct frame *f, size_t ne)
{
	params_reinit(&rgrad->params, f, ne);
}

static void rgrad_deinit(struct recv_fit_rgrad *rgrad)
{
	params_deinit(&rgrad->params);
}

static void rgrad_set(struct recv_fit_rgrad *rgrad, size_t dim, size_t ne,
		      const struct recv_fit_kkt *kkt,
		      const struct recv_fit_resid *resid)
{
	assert(!kkt->factored);

	size_t n = dim + ne;
	enum blas_uplo uplo = MLOGIT_COV_UPLO == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER;

	if (!n)
		return;

	blas_dsymv(uplo, n, 2.0, kkt->matrix, n, resid->params.all, 1, 0.0, rgrad->params.all, 1);
}

static enum recv_fit_task primal_dual_step(struct recv_fit *fit)
{
	size_t dim = recv_model_dim(&fit->model);
	size_t ne = recv_fit_constr_count(fit);
	size_t n = dim + ne;

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

		resid_set(&fit->cur->resid, &fit->cur->loglik, &fit->constr, &fit->cur->params);
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
	kkt_set(&fit->kkt, &fit->cur->loglik, &fit->constr);
	rgrad_set(&fit->rgrad, dim, ne, &fit->kkt, &fit->cur->resid); // TODO: remove this?

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
		   const struct messages *xmsgs,
		   const struct messages *ymsgs,
		   const struct recv_fit_ctrl *ctrl)
{
	assert(fit);
	assert(xmsgs);
	assert(f);
	assert(messages_max_to(xmsgs) < frame_send_count(f));
	assert(messages_max_from(xmsgs) < frame_recv_count(f));
	assert(!ymsgs || (messages_max_to(ymsgs) < frame_send_count(f)));
	assert(!ymsgs || (messages_max_from(ymsgs) < frame_recv_count(f)));
	assert(!ctrl || recv_fit_ctrl_valid(ctrl));

	struct history *h = frame_history(f);
	size_t ne = 0;

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
	constr_init(&fit->constr);
	eval_init(&fit->eval[0], &fit->model, ne);
	eval_init(&fit->eval[1], &fit->model, ne);
	fit->prev = &fit->eval[0];
	fit->cur = &fit->eval[1];
	kkt_init(&fit->kkt, f, ne);
	search_init(&fit->search, f, ne);
	rgrad_init(&fit->rgrad, f, ne);

	// evaluate the model at zero
	eval_set(fit->cur, NULL, fit->xmsgs, fit->ymsgs, fit->frame, &fit->model, 2);

	fit->dev0 = recv_fit_dev(fit);
}


enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct recv_fit_params *params0)
{
	if (params0 != NULL) {
		eval_set(fit->cur, params0, fit->xmsgs, fit->ymsgs, fit->frame,
			 &fit->model, 2);
	}

	resid_set(&fit->cur->resid, &fit->cur->loglik, &fit->constr, &fit->cur->params);
	kkt_set(&fit->kkt, &fit->cur->loglik, &fit->constr);
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
	constr_deinit(&fit->constr);
	recv_model_deinit(&fit->model);
}

size_t recv_fit_constr_count(const struct recv_fit *fit)
{
	return constr_count(&fit->constr);
}

void recv_fit_get_constr(const struct recv_fit *fit, size_t i,
			 const double **weightsp, double *valp)
{
	size_t dim = recv_model_dim(&fit->model);
	constr_get(&fit->constr, dim, i, weightsp, valp);
}

const double *recv_fit_constrs(const struct recv_fit *fit)
{
	return fit->constr.wts;
}

const double *recv_fit_constr_vals(const struct recv_fit *fit)
{
	return fit->constr.vals;
}

static void _recv_fit_add_constrs(struct recv_fit *fit)
{
	const struct frame *f = fit->frame;
	size_t ne = recv_fit_constr_count(fit);

	eval_reinit(&fit->eval[0], f, ne);
	eval_reinit(&fit->eval[1], f, ne);
	kkt_reinit(&fit->kkt, f, ne);
	search_reinit(&fit->search, f, ne);
	rgrad_reinit(&fit->rgrad, f, ne);
}

void recv_fit_add_constr(struct recv_fit *fit, const double *weights,
			 double val)
{
	size_t dim = recv_model_dim(&fit->model);
	constr_add(&fit->constr, dim, weights, val);
	_recv_fit_add_constrs(fit);
}

void recv_fit_add_constr_set(struct recv_fit *fit, size_t i, double val)
{
	assert(i < recv_model_dim(&fit->model));
	
	size_t dim = recv_model_dim(&fit->model);
	constr_add_set(&fit->constr, dim, i, val);
	_recv_fit_add_constrs(fit);
}

void recv_fit_add_constr_eq(struct recv_fit *fit, size_t i1, size_t i2)
{
	assert(i1 < recv_model_dim(&fit->model));
	assert(i2 < recv_model_dim(&fit->model));
	assert(i1 != i2);

	size_t dim = recv_model_dim(&fit->model);
	constr_add_eq(&fit->constr, dim, i1, i2);
	_recv_fit_add_constrs(fit);
}


size_t recv_fit_add_constr_identify(struct recv_fit *fit)
{
	size_t iadd, nadd = 0;
	const struct recv_loglik *ll = &fit->cur->loglik;
	const struct recv_model *m = recv_loglik_model(ll);
	size_t i, n = recv_model_dim(m);

	if (n == 0)
		goto out;

	// copy the information matrix at 0
	double *imat = xmalloc(n * n * sizeof(imat[0]));
	double *imat_buf = xmalloc(n * (n + 1) / 2 * sizeof(imat_buf[0]));
	double *scale = xmalloc(n * sizeof(scale[0]));
	memset(imat, 0, n * n * sizeof(imat[0]));
	memset(imat_buf, 0, n * (n + 1) / 2 * sizeof(imat_buf[0]));

	double df = (double)recv_loglik_count(ll);
	recv_loglik_axpy_imat(1.0/df, ll, imat_buf);

	double *dst = imat;
	const double *src = imat_buf;

	if (MLOGIT_COV_UPLO == BLAS_UPPER) {
		for (i = 0; i < n; i++) {
			size_t rowlen = n - i;

			memcpy(dst + i, src, rowlen * sizeof(*dst));
			scale[i] = 1.0 / (dst[i] < fit->ctrl.vartol ? 1.0 : sqrt(dst[i]));

			dst += n;
			src += rowlen;
		}
	} else {
		for (i = 0; i < n; i++) {
			size_t rowlen = i + 1;

			memcpy(dst, src, rowlen * sizeof(*dst));
			scale[i] = 1.0 / (dst[i] < fit->ctrl.vartol ? 1.0 : sqrt(dst[i]));

			dst += n;
			src += rowlen;
		}
	}

	// scale the information matrix
	for (i = 0; i < n; i++) {
		blas_dscal(n, scale[i], imat + i * n, 1); // columns
	}
	for (i = 0; i < n; i++) {
		blas_dscal(n, scale[i], imat + i, n); // rows
	}


	// compute the eigendecomposition of the information matrix
	enum lapack_eigjob jobz = LA_EIG_VEC;
	enum blas_uplo uplo = MLOGIT_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER;
	size_t lwork, liwork;
	lwork = lapack_dsyevd_lwork(jobz, n, &liwork);
	double *work = xmalloc(lwork * sizeof(work[0]));
	ptrdiff_t *iwork = xmalloc(liwork * sizeof(iwork[0]));
	double *evals = xmalloc(n * sizeof(evals[0]));
	ptrdiff_t info = lapack_dsyevd(jobz, uplo, n, imat, n, evals, work, lwork, iwork, liwork);
	assert(info == 0);

	// change back to the original scale
	for (i = 0; i < n; i++) {
		blas_dscal(n, scale[i], imat + i, n);
	}

	// compute the dimension of the nullspace
	size_t nulldim;
	for (nulldim = 0; nulldim < n; nulldim++) {
		if (evals[nulldim] > fit->ctrl.eigtol)
			break;
	}
	// at this point, the first nulldim columns of imat store its nullspace

	if (nulldim == 0)
		goto cleanup_eig_imat;


	// compute the projection of the constraints onto the nullspace
	const struct recv_fit_constr *ce = &fit->constr;
	size_t ne = constr_count(ce);

	if (ne > 0) {
		double *proj = xmalloc(nulldim * ne * sizeof(proj[0]));
		blas_dgemm(BLAS_TRANS, BLAS_NOTRANS, nulldim, ne, n, 1.0, imat, n, ce->wts, n, 0.0, proj, nulldim);

		// compute the nullspace of the projection
		double *s = xmalloc(MIN(ne, nulldim) * sizeof(s[0]));
		double *u = xmalloc(nulldim * nulldim * sizeof(u[0]));
		double *vt = xmalloc(ne * ne * sizeof(vt[0]));

		enum lapack_svdjob svd_jobz = LA_SVD_ALL;
		lwork = lapack_dgesdd_lwork(svd_jobz, nulldim, ne, &liwork);
		work = xrealloc(work, lwork * sizeof(work[0]));
		iwork = xrealloc(iwork, liwork * sizeof(iwork[0]));
		ptrdiff_t svd_info = lapack_dgesdd(svd_jobz, nulldim, ne, proj, nulldim, s, u, nulldim, vt, ne, work, lwork, iwork);
		assert(svd_info == 0);

		// compute the rank
		size_t rank, maxrank = MIN(ne, nulldim);
		for (rank = maxrank; rank > 0; rank--) {
			if (s[rank - 1] > fit->ctrl.svtol * fabs(s[maxrank - 1]))
				break;
		}

		// at this point the last nulldim - rank columns of u hold the coefficients of
		// the orthogonal complement to the constraint matrix in the nullspace


		// compute a basis for the orthogonal complement; these are the consstraints
		nadd = nulldim - rank;
		double *compl = xmalloc(n * nadd * sizeof(compl[0]));
		blas_dgemm(BLAS_NOTRANS, BLAS_NOTRANS, n, nadd, nulldim, 1.0, imat, n, u + rank, nulldim, 0.0, compl, n);

		for (iadd = 0; iadd < nadd; iadd++) {
			constr_add(&fit->constr, n, compl + iadd * n, 0.0);
		}

		free(compl);
		free(vt);
		free(u);
		free(s);
		free(proj);
	} else {
		// add the nullspace of the information matrix as constraints
		nadd = nulldim;
		for (iadd = 0; iadd < nadd; iadd++) {
			constr_add(&fit->constr, n, imat + iadd * n, 0.0);
		}
	}
	_recv_fit_add_constrs(fit);

cleanup_eig_imat:
	free(evals);
	free(iwork);
	free(work);
	free(scale);
	free(imat_buf);
	free(imat);
out:
	return nadd;
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

size_t recv_fit_coefs_count(const struct recv_fit *fit)
{
	const struct recv_model *m = &fit->model;
	return recv_model_dim(m);
}
