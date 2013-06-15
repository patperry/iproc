#include "port.h"
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "matrixutil.h"
#include "lapack.h"
#include "xalloc.h"
#include "newton.h"


static void params_init(struct newton_params *p, size_t dim, size_t nc);
static void params_deinit(struct newton_params *p);
static double *params_all(const struct newton_params *p);
static double *params_pri(const struct newton_params *p);
static double *params_dual(const struct newton_params *p);
static void params_set(struct newton_params *p, const double *pri,
		       const double *dual, size_t dim, size_t nc);
static void params_assign(struct newton_params *dst,
			  const struct newton_params *src, size_t dim,
			  size_t nc);
static void params_axpy(double alpha, const struct newton_params *x,
			struct newton_params *y, size_t dim, size_t nc);

static void eval_init(struct newton_eval *e, size_t dim, size_t nc);
static void eval_deinit(struct newton_eval *e);

static void kkt_init(struct newton_kkt *kkt, size_t dim, size_t nc);
static void kkt_deinit(struct newton_kkt *kkt);
static void kkt_set(struct newton_kkt *kkt, const double *hess, size_t dim,
		    enum blas_uplo uplo, const struct constr *c);
static int kkt_solve(struct newton_kkt *kkt, struct newton_params *x,
		     size_t dim, size_t nc);



static void params_init(struct newton_params *p, size_t dim, size_t nc)
{
	p->pri = xcalloc(dim + nc, sizeof(double));
	p->dual = p->pri + dim;
}


static void params_deinit(struct newton_params *p)
{
	free(p->pri);
}


static double *params_all(const struct newton_params *p)
{
	return p->pri;
}


static double *params_pri(const struct newton_params *p)
{
	return p->pri;
}


static double *params_dual(const struct newton_params *p)
{
	return p->dual;
}

static void params_set(struct newton_params *p, const double *pri,
		       const double *dual, size_t dim, size_t nc)
{
	memset(params_all(p), 0, (dim + nc) * sizeof(double));
	if (pri)
		memcpy(params_pri(p), pri, dim * sizeof(double));
	if (dual)
		memcpy(params_dual(p), dual, nc * sizeof(double));
}

static void params_assign(struct newton_params *dst,
			  const struct newton_params *src, size_t dim,
			  size_t nc)
{
	if (src) {
		params_set(dst, src->pri, src->dual, dim, nc);
	} else {
		params_set(dst, NULL, NULL, dim, nc);
	}
}

static void params_axpy(double alpha, const struct newton_params *x,
			struct newton_params *y, size_t dim, size_t nc)
{
	size_t n = dim + nc;
	blas_daxpy(n, alpha, params_all(x), 1, params_all(y), 1);
}

static void eval_init(struct newton_eval *e, size_t dim, size_t nc)
{
	params_init(&e->params, dim, nc);
	e->val = NAN;
	e->grad = xmalloc(dim * sizeof(double));

	params_init(&e->resid, dim, nc);
	e->resid_norm = NAN;
	e->feasible = 0;
	e->indomain = 0;
}


static void eval_set_resid(struct newton_eval *e, size_t dim,
			   const struct constr *c)
{
	size_t nc = c ? constr_count(c) : 0;
	size_t n = dim + nc;
	const double *A = c ? constr_all_wts(c) : NULL;
	const double *b = c ? constr_all_vals(c) : NULL;
	const double *x = params_pri(&e->params);
	const double *nu = params_dual(&e->params);
	const double *grad = e->grad;
	double *rdual = params_pri(&e->resid);
	double *rpri = params_dual(&e->resid);
	double *r = params_all(&e->resid);
	double norm;

	if (!dim)
		goto set_norm;

	/* dual residual: grad(f) + A * nu */
	memcpy(rdual, grad, dim * sizeof(double));
	blas_dgemv(BLAS_NOTRANS, dim, nc, 1.0, A, dim, nu, 1, 1.0, rdual, 1);

	/* primal residual: A' * x - b */
	memcpy(rpri, b, nc * sizeof(double));
	blas_dgemv(BLAS_TRANS, dim, nc, 1.0, A, dim, x, 1, -1.0, rpri, 1);

set_norm:
	/* compute ||r|| */
	norm = blas_dnrm2(n, r, 1);
	e->resid_norm = norm;
}


static void eval_set(struct newton_eval *e, double val, const double *grad,
		     size_t dim, const struct constr *c)
{
	size_t nc = c ? constr_count(c) : 0;

	e->val = val;
	memcpy(e->grad, grad, dim * sizeof(double));

	/* TODO: better feasible test, check for nans in resid, grad */
	eval_set_resid(e, dim, c);
	e->feasible = blas_dnrm2(nc, params_dual(&e->resid), 1) <= nc * 1e-8;
	e->indomain = isfinite(val) && isfinite(e->resid_norm);
}



static void eval_deinit(struct newton_eval *e)
{
	params_deinit(&e->resid);
	free(e->grad);
	params_deinit(&e->params);
}


static void kkt_init(struct newton_kkt *kkt, size_t dim, size_t nc)
{
	size_t n = dim + nc;
	size_t lwork = lapack_dsysv_lwork(n);

	kkt->factored = 0;
	kkt->matrix = xmalloc(n * n * sizeof(double));
	kkt->ldl_ipiv = xmalloc(n * sizeof(ptrdiff_t));
	kkt->ldl_work = xmalloc(lwork * sizeof(double));
	kkt->ldl_lwork = lwork;
}


static void kkt_deinit(struct newton_kkt *kkt)
{
	free(kkt->ldl_work);
	free(kkt->ldl_ipiv);
	free(kkt->matrix);
}


static void kkt_set(struct newton_kkt *kkt, const double *hess, size_t dim,
		    enum blas_uplo uplo, const struct constr *c)
{
	size_t nc = c ? constr_count(c) : 0;
	size_t n = dim + nc;
	enum blas_uplo f77uplo = (uplo == BLAS_UPPER) ? BLAS_LOWER : BLAS_UPPER;

	kkt->factored = 0;

	memset(kkt->matrix, 0, n * n * sizeof(double));
	kkt->uplo = uplo;

	/* k11 */
	packed_dsctr(f77uplo, dim, hess, kkt->matrix, n);

	if (!nc)
		return;

	if (uplo == BLAS_UPPER) { /* k12 */
		matrix_dtrans(dim, nc, c->wts, dim, kkt->matrix + dim, n);
	} else { /* k21 */
		lapack_dlacpy(LA_COPY_ALL, dim, nc, c->wts, dim,
			      kkt->matrix + dim * n, n);
	}

	kkt->uplo = uplo;
}


static int kkt_solve(struct newton_kkt *kkt, struct newton_params *x,
		     size_t dim, size_t nc)
{
	assert(!kkt->factored);

	size_t n = dim + nc;
	enum blas_uplo f77uplo = (kkt->uplo == BLAS_UPPER) ? BLAS_LOWER : BLAS_UPPER;
	ptrdiff_t info = 0;
	int err = 0;

	if (n) {
		info = lapack_dsysv(f77uplo, n, 1, kkt->matrix, n,
				    kkt->ldl_ipiv, params_all(x), n,
				    kkt->ldl_work, kkt->ldl_lwork);
	}

	kkt->factored = 1;

	if (info) {
		err = EDOM;
		errno = err;
	}

	return -err;
}


static int search_set(struct newton_params *search,
		      const struct newton_params *resid, struct newton_kkt *kkt,
		      size_t dim, size_t nc)
{
	size_t n = dim + nc;
	double *s = params_all(search);
	int err;

	memcpy(s, params_all(resid), n * sizeof(double));
	blas_dscal(n, -1.0, s, 1);
	err = kkt_solve(kkt, search, dim, nc);

	return err;
}


void newton_init(struct newton *opt, size_t dim,
		 const struct constr *constr,
		 const struct newton_ctrl *ctrl)
{
	assert(!constr || constr_dim(constr) == dim);
	assert(!ctrl || newton_ctrl_valid(ctrl));

	size_t nc = constr ? constr_count(constr) : 0;

	if (!ctrl) {
		opt->ctrl = NEWTON_CTRL0;
	} else {
		opt->ctrl = *ctrl;
	}
	opt->dim = dim;
	opt->constr = constr;

	eval_init(&opt->eval[0], dim, nc);
	eval_init(&opt->eval[1], dim, nc);
	opt->cur = NULL;
	opt->next = NULL;

	kkt_init(&opt->kkt, dim, nc);
	params_init(&opt->search, dim, nc);
	opt->step = NAN;
}


void newton_deinit(struct newton *opt)
{
	params_deinit(&opt->search);
	kkt_deinit(&opt->kkt);
	eval_deinit(&opt->eval[1]);
	eval_deinit(&opt->eval[0]);
}


enum newton_task newton_start(struct newton *opt, const double *x0,
			      double f0, const double *grad0,
			      const double *duals)
{
	size_t dim = newton_dim(opt);
	const struct constr *c = newton_constr(opt);
	size_t nc = c ? constr_count(c) : 0;
	const struct newton_ctrl *ctrl = newton_ctrl(opt);
	enum newton_task task;

	opt->cur = &opt->eval[0];
	opt->next = NULL;

	params_set(&opt->cur->params, x0, duals, dim, nc);
	eval_set(opt->cur, f0, grad0, dim, c);

	/* ensure x0 is in domain */
	if (!opt->cur->indomain) {
		task = NEWTON_ERR_DOM;

	/* stop early if residual is below tolerance */
	} else if (opt->cur->feasible && opt->cur->resid_norm <= ctrl->gtol) {
		task = NEWTON_CONV;

	/* compute first newton step */
	} else {
		task = NEWTON_HESS;
	}

	return task;
}


const double *newton_next(const struct newton *opt)
{
	assert(opt->next);
	return params_pri(&opt->next->params);
}


enum newton_task newton_step(struct newton *opt, double f, const double *grad)
{
	assert(opt->next);

	size_t dim = newton_dim(opt);
	const struct constr *c = newton_constr(opt);
	size_t nc = c ? constr_count(c) : 0;
	const struct newton_ctrl *ctrl = newton_ctrl(opt);
	struct newton_eval *cur = opt->cur;
	struct newton_eval *next = opt->next;
	double beta = 0.5;
	double alpha = 0.1;
	double t = opt->step;
	double t1 = beta * t;
	enum newton_task task;

	opt->ls_iter++;

	if (isnan(f) || !grad)
		goto shrink;

	eval_set(next, f, grad, dim, c);

	if (next->indomain) {
		if (next->feasible && next->resid_norm <= ctrl->gtol) {
			task = NEWTON_CONV;
			goto out;
		} else if (next->resid_norm <= (1 - alpha * t) * cur->resid_norm) {
			task = NEWTON_HESS;
			goto out;
		}
	}

shrink:
	if (t1 <= ctrl->xtol) {
		task = NEWTON_ERR_XTOL;
	} else if (opt->ls_iter == ctrl->ls_maxit) {
		task = NEWTON_ERR_LNSRCH;
	} else {
		opt->step = t1;
		params_assign(&next->params, &cur->params, dim, nc);
		params_axpy(opt->step, &opt->search, &next->params, dim, nc);
		task = NEWTON_STEP;
	}

out:
	if (task == NEWTON_CONV || task == NEWTON_HESS) {
		opt->cur = next;
	}
	if (task != NEWTON_STEP) {
		opt->next = NULL;
	}
	return task;
}


enum newton_task newton_set_hess(struct newton *opt, const double *hess,
				 enum blas_uplo uplo)
{
	assert(!opt->next);

	size_t dim = newton_dim(opt);
	const struct constr *c = newton_constr(opt);
	size_t nc = c ? constr_count(c) : 0;
	int err;

	/* compute search direction */
	kkt_set(&opt->kkt, hess, dim, uplo, c);
	err = search_set(&opt->search, &opt->cur->resid, &opt->kkt, dim, nc);
	if (err == -EDOM)
		return NEWTON_ERR_HESS;
	assert(err == 0);

	/* start a new line search */
	opt->step = 1.0;
	opt->ls_iter = 0;

	/* compute next x */
	opt->next = opt->cur == &opt->eval[0] ? &opt->eval[1] : &opt->eval[0];
	params_assign(&opt->next->params, &opt->cur->params, dim, nc);
	params_axpy(opt->step, &opt->search, &opt->next->params, dim, nc);

	return NEWTON_STEP;
}

const char *newton_errmsg(enum newton_task task)
{
	switch (task) {
	case NEWTON_STEP:
	case NEWTON_HESS:
		return "optimization in progress";
	case NEWTON_ERR_LNSRCH:
		return "linesearch failed";
	case NEWTON_ERR_XTOL:
		return "step size is less than tolerance";
	case NEWTON_ERR_HESS:
		return "Hessian is singular";
	case NEWTON_ERR_DOM:
		return "initial point is not in domain";
	case NEWTON_CONV:
		return NULL;
	}

	assert(0);
	return NULL;
}

/* current values */
const double *newton_params(const struct newton *opt)
{
	assert(opt->cur);
	return params_pri(&opt->cur->params);
}

const double *newton_duals(const struct newton *opt)
{
	assert(opt->cur);
	return params_dual(&opt->cur->params);
}

double newton_val(const struct newton *opt)
{
	assert(opt->cur);
	return opt->cur->val;
}

const double *newton_grad(const struct newton *opt)
{
	assert(opt->cur);
	return opt->cur->grad;
}

double newton_grad_norm(const struct newton *opt)
{
	assert(opt->cur);
	size_t n = opt->dim;
	const double *g = params_pri(&opt->cur->resid);
	double nrm = blas_dnrm2(n, g, 1);
	return nrm;
}

double newton_resid_norm(const struct newton *opt)
{
	assert(opt->cur);
	return opt->cur->resid_norm;
}

double newton_step_size(const struct newton *opt)
{
	assert(opt->cur);
	return opt->step;
}

const double *newton_search(const struct newton *opt)
{
	assert(opt->cur);
	return params_pri(&opt->search);
}

int newton_ctrl_valid(const struct newton_ctrl *ctrl)
{
	assert(ctrl);

	if (!(ctrl->gtol > 0)) {
		return 0;
	} else {
		return linesearch_ctrl_valid(&ctrl->ls);
	}
}

