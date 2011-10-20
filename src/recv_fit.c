#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sblas.h"
#include "xalloc.h"
#include "linesearch.h"
#include "recv_fit.h"

static void constr_init(struct recv_fit_constr *c)
{
	c->wts = NULL;
	c->wt_inds = NULL;
	c->wt_cols = NULL;
	c->vals = NULL;
	c->names = NULL;
	c->n = 0;
	c->nmax = 0;
	c->wt_nzmax = 0;
}

static void constr_deinit(struct recv_fit_constr *c)
{
	size_t i, n = c->n;

	for (i = 0; i < n; i++) {
		free(c->names[i]);
	}
	free(c->names);
	free(c->vals);
	free(c->wt_cols);
	free(c->wt_inds);
	free(c->wts);
}

static void constr_grow(struct recv_fit_constr *c)
{
	if (c->n == c->nmax) {
		size_t nmax = ARRAY_GROW(c->nmax, SIZE_MAX);
		c->wt_cols = xrealloc(c->wt_cols,
				      (nmax + 1) * sizeof(c->wt_cols[0]));
		c->vals = xrealloc(c->vals, nmax * sizeof(c->vals[0]));
		c->names = xrealloc(c->names, nmax * sizeof(c->names[0]));
		c->nmax = nmax;

		if (!c->n)
			c->wt_cols[0] = 0;
	}
}

static void constr_grow_nz(struct recv_fit_constr *c, size_t delta)
{
	size_t nz0 = c->n ? c->wt_cols[c->n] : 0;
	assert(delta <= SIZE_MAX - nz0);

	size_t nz = nz0 + delta;
	size_t nzmax = c->wt_nzmax;

	if (nz <= nzmax)
		return;

	while (nz > nzmax) {
		nzmax = ARRAY_GROW(nzmax, SIZE_MAX);
	}

	c->wts = xrealloc(c->wts, nzmax * sizeof(c->wts[0]));
	c->wt_inds = xrealloc(c->wt_inds, nzmax * sizeof(c->wt_inds[0]));
	c->wt_nzmax = nzmax;
}

static size_t constr_count(const struct recv_fit_constr *c)
{
	return c->n;
}

static void constr_add(struct recv_fit_constr *c, const double *wts,
		       const size_t *ind, size_t nz,
		       double val, const char *name)
{
	constr_grow(c);
	constr_grow_nz(c, nz);

	size_t n = c->n;
	memcpy(c->wts + c->wt_cols[n], wts, nz * sizeof(wts[0]));
	memcpy(c->wt_inds + c->wt_cols[n], ind, nz * sizeof(ind[0]));
	c->wt_cols[n+1] = c->wt_cols[n] + nz;
	c->vals[n] = val;
	c->names[n] = xstrdup(name);
	c->n = n + 1;
}

static void constr_add_set(struct recv_fit_constr *c, size_t dim, size_t nc,
			   size_t i, size_t ic, double val)
{
	assert(c);
	assert(i < dim);
	assert(ic < nc);
	assert(isfinite(val));

	const char *fmt = "Set(%zd,%zd,%g)";
	char buf;
	size_t len = snprintf(&buf, 1, fmt, i + 1, ic + 1, val);
	char *name = xmalloc((len + 1) * sizeof(char));
	snprintf(name, len + 1, fmt, i + 1, ic + 1, val);

	constr_grow(c);
	constr_grow_nz(c, 1);

	size_t n = c->n;
	c->wts[c->wt_cols[n]] = 1.0;
	c->wt_inds[c->wt_cols[n]] = i + ic * dim;
	c->wt_cols[n+1] = c->wt_cols[n] + 1;
	c->vals[n] = val;
	c->names[n] = name;
	c->n = n + 1;
}

static void constr_add_eq(struct recv_fit_constr *c, size_t dim, size_t nc,
			  size_t i1, size_t ic1, size_t i2, size_t ic2)
{
	assert(c);
	assert(i1 < dim);
	assert(i2 < dim);
	assert(ic1 < nc);
	assert(ic2 < nc);
	assert(i1 != i2 || ic1 != ic2);

	const char *fmt = "Eq((%zd,%zd),(%zd,%zd))";
	char buf;
	size_t len = snprintf(&buf, 1, fmt, i1 + 1, ic1 + 1, i2 + 1, ic2 + 1);
	char *name = xmalloc((len + 1) * sizeof(char));
	snprintf(name, len + 1, fmt, i1 + 1, ic1 + 1, i2 + 1, ic2 + 1);

	constr_grow(c);
	constr_grow_nz(c, 2);

	size_t n = c->n;
	c->wts[c->wt_cols[n]] = +1.0;
	c->wt_inds[c->wt_cols[n]] = i1 + ic1 * dim;
	c->wts[c->wt_cols[n] + 1] = -1.0;
	c->wt_inds[c->wt_cols[n]] = i2 + ic2 * dim;
	c->wt_cols[n+1] = c->wt_cols[n] + 2;
	c->vals[n] = 0.0;
	c->names[n] = name;
	c->n = n + 1;
}

static void constr_get(const struct recv_fit_constr *c, size_t i,
		       const double **weightsp, const size_t **indp,
		       size_t *nzp, double *valp, const char **namep)
{
	assert(i < constr_count(c));

	const double *wts = c->wts + c->wt_cols[i];
	const size_t *ind = c->wt_inds + c->wt_cols[i];
	size_t nz = c->wt_cols[i+1] - c->wt_cols[i];
	double val = c->vals[i];
	const char *name = c->names[i];

	if (weightsp)
		*weightsp = wts;
	if (indp)
		*indp = ind;
	if (nzp)
		*nzp = nz;
	if (valp)
		*valp = val;
	if (namep)
		*namep = name;

}

static void resid_init(struct recv_fit_resid *resid, ssize_t nc, ssize_t dim,
		       ssize_t ne)
{
	vector_init(&resid->vector, dim * nc + ne);
}

static void resid_deinit(struct recv_fit_resid *resid)
{
	vector_deinit(&resid->vector);
}

static void resid_add_constrs(struct recv_fit_resid *resid, ssize_t n)
{
	assert(n >= 0);

	ssize_t m0 = vector_dim(&resid->vector);
	ssize_t m = m0 + n;
	vector_reinit(&resid->vector, m);
}

static void resid_set(struct recv_fit_resid *resid,
		      const struct recv_loglik *ll,
		      const struct recv_fit_constr *ce,
		      const struct vector *params)
{
	struct vector *r = &resid->vector;
	const struct recv_model *m = recv_loglik_model(ll);
	size_t ntot = recv_loglik_count_sum(ll);
	size_t dim = recv_model_dim(m);
	size_t ic, nc = recv_model_cohort_count(m);
	size_t nce = constr_count(ce);

	/* r1 is the dual residual: grad(f) + ce * nu,
	 *   where grad(f) = grad(nll)
	 *                 = -(score)
	 */
	struct vector r1 = vector_slice(r, 0, dim * nc);
	struct vector duals = vector_slice(params, dim * nc, nce);

	vector_fill(&r1, 0.0);
	for (ic = 0; ic < nc; ic++) {
		const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
		const struct vector *score = &info->score;
		size_t n = recv_loglik_count(ll, ic);

		struct vector r1c = vector_slice(&r1, ic * dim, dim);
		vector_axpy(-((double)n) / ntot, score, &r1c);
	}

	sblas_dcscmv(BLAS_NOTRANS, dim * nc, nce, 1.0, ce->wts, ce->wt_inds,
		     ce->wt_cols, vector_to_ptr(&duals), 1.0,
		     vector_to_ptr(&r1));

	// r2 is the primal residual: ce' * x - be
	struct vector r2 = vector_slice(r, dim * nc, nce);
	struct vector primals = vector_slice(params, 0, dim * nc);

	blas_dcopy(nce, ce->vals, 1, vector_to_ptr(&r2), 1);
	sblas_dcscmv(BLAS_TRANS, dim * nc, nce, 1.0, ce->wts, ce->wt_inds,
		     ce->wt_cols, vector_to_ptr(&primals), -1.0,
		     vector_to_ptr(&r2));

	// compute ||r||^2
	resid->norm2 = vector_norm2(r);
}

static void _eval_setup(struct recv_fit_eval *eval, ssize_t dim, ssize_t nc,
			ssize_t ne)
{
	struct vector coefs_data = vector_slice(&eval->params, 0, dim * nc);
	eval->coefs = matrix_make(&coefs_data, dim, nc);
	eval->duals = vector_slice(&eval->params, dim * nc, ne);
}

static void eval_init(struct recv_fit_eval *eval, struct recv_model *m,
		      ssize_t ne)
{
	ssize_t dim = recv_model_dim(m);
	ssize_t nc = recv_model_cohort_count(m);

	vector_init(&eval->params, nc * dim + ne);
	_eval_setup(eval, dim, nc, ne);

	recv_loglik_init(&eval->loglik, m);
	resid_init(&eval->resid, nc, dim, ne);
}

static void eval_deinit(struct recv_fit_eval *eval)
{
	resid_deinit(&eval->resid);
	recv_loglik_deinit(&eval->loglik);
	vector_deinit(&eval->params);
}

static void eval_add_constrs(struct recv_fit_eval *eval, ssize_t n)
{
	assert(n >= 0);

	ssize_t dim = matrix_nrow(&eval->coefs);
	ssize_t nc = matrix_ncol(&eval->coefs);
	ssize_t ne0 = vector_dim(&eval->duals);
	ssize_t ne = ne0 + n;

	vector_reinit(&eval->params, nc * dim + ne);
	_eval_setup(eval, dim, nc, ne);
	resid_add_constrs(&eval->resid, n);
}

static void _eval_update(struct recv_fit_eval *eval,
			 const struct recv_fit_constr *ce,
			 const struct messages *xmsgs,
			 const struct messages *ymsgs,
			 struct frame *frame, struct recv_model *model)
{
	// clear frame; set model coefs
	frame_clear(frame);
	recv_model_set_coefs(model, &eval->coefs);

	// set loglik
	recv_loglik_clear(&eval->loglik);
	struct messages_iter xit = messages_iter_make(xmsgs);
	struct messages_iter yit = messages_iter_make(ymsgs);
	double xt = -INFINITY;
	double yt = -INFINITY;
	const struct message *msg;
	ssize_t i, n;

	if (messages_iter_advance(&xit)) {
		xt = MESSAGES_TIME(xit);
	} else {
		xt = INFINITY;
	}

	while (messages_iter_advance(&yit)) {
		yt = MESSAGES_TIME(yit);

		while (xt < yt) {
			frame_advance(frame, xt);

			n = MESSAGES_COUNT(xit);
			for (i = 0; i < n; i++) {
				msg = MESSAGES_VAL(xit, i);
				frame_add(frame, msg);
			}

			if (messages_iter_advance(&xit)) {
				xt = MESSAGES_TIME(xit);
			} else {
				xt = INFINITY;
			}
		}

		frame_advance(frame, yt);
		n = MESSAGES_COUNT(yit);
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(yit, i);
			recv_loglik_add(&eval->loglik, frame, msg);
		}
	}

	// set resid
	resid_set(&eval->resid, &eval->loglik, ce, &eval->params);

	// set in_domain
	eval->in_domain = isfinite(eval->resid.norm2);

}

static void eval_set(struct recv_fit_eval *eval,
		     const struct recv_fit_constr *ce,
		     const struct matrix *coefs, const struct vector *duals,
		     const struct messages *xmsgs,
		     const struct messages *ymsgs,
		     struct frame *frame, struct recv_model *model)
{
	size_t dim = recv_model_dim(model);
	size_t nc = recv_model_cohort_count(model);
	size_t nce = constr_count(ce);

	struct vector pcoefs_data = vector_slice(&eval->params, 0, dim * nc);
	struct matrix pcoefs = matrix_make(&pcoefs_data, dim, nc);

	if (coefs) {
		matrix_assign_copy(&pcoefs, BLAS_NOTRANS, coefs);
	} else {
		matrix_fill(&pcoefs, 0.0);
	}

	struct vector pduals = vector_slice(&eval->params, dim * nc, nce);

	if (duals) {
		vector_assign_copy(&pduals, duals);
	} else {
		vector_fill(&pduals, 0.0);
	}

	_eval_update(eval, ce, xmsgs, ymsgs, frame, model);
}

static void eval_step(struct recv_fit_eval *eval,
		      const struct recv_fit_constr *ce,
		      const struct vector *params0, double scale,
		      const struct vector *dir, const struct messages *xmsgs,
		      const struct messages *ymsgs,
		      struct frame *frame, struct recv_model *model)
{
	vector_assign_copy(&eval->params, params0);
	vector_axpy(scale, dir, &eval->params);
	_eval_update(eval, ce, xmsgs, ymsgs, frame, model);
}

static void kkt_init(struct recv_fit_kkt *kkt, ssize_t nc, ssize_t dim,
		     ssize_t nce)
{
	ssize_t n = nc * dim + nce;

	kkt->factored = false;
	kkt->uplo = BLAS_LOWER;
	matrix_init(&kkt->matrix, n, n);
	ldlfac_init(&kkt->ldl, n);
}

static void kkt_add_constrs(struct recv_fit_kkt *kkt, ssize_t n)
{
	assert(n >= 0);
	ssize_t m = matrix_nrow(&kkt->matrix) + n;

	kkt->factored = false;
	matrix_reinit(&kkt->matrix, m, m);
	ldlfac_reinit(&kkt->ldl, m);
}

static void kkt_deinit(struct recv_fit_kkt *kkt)
{
	ldlfac_deinit(&kkt->ldl);
	matrix_deinit(&kkt->matrix);
}

static void kkt_set(struct recv_fit_kkt *kkt, const struct recv_loglik *ll,
		    const struct recv_fit_constr *ce)
{
	struct matrix *k = &kkt->matrix;
	const struct recv_model *m = recv_loglik_model(ll);
	size_t dim = recv_model_dim(m);
	size_t ic, nc = recv_model_cohort_count(m);
	size_t ntot = recv_loglik_count_sum(ll);
	size_t nce = constr_count(ce);

	matrix_fill(k, 0.0);

	// k11
	for (ic = 0; ic < nc; ic++) {
		struct recv_loglik_info *info = recv_loglik_info(ll, ic);
		ssize_t n = recv_loglik_count(ll, ic);
		struct matrix h = matrix_slice(k, ic * dim, ic * dim, dim, dim);
		matrix_assign_copy(&h, BLAS_NOTRANS, &info->imat);
		matrix_scale(&h, ((double)n) / ntot);
	}

	// k12, k21
	struct matrix k12 = matrix_slice(k, 0, nc * dim, nc * dim, nce);
	struct matrix k21 = matrix_slice(k, nc * dim, 0, nce, nc * dim);

	matrix_fill(&k12, 0.0);
	sblas_dcscsctr(BLAS_NOTRANS, ce->n, ce->wts, ce->wt_inds, ce->wt_cols,
		       matrix_to_ptr(&k12), matrix_lda(&k12));

	matrix_fill(&k21, 0.0);
	sblas_dcscsctr(BLAS_TRANS, ce->n, ce->wts, ce->wt_inds, ce->wt_cols,
		       matrix_to_ptr(&k21), matrix_lda(&k21));

	kkt->factored = false;
}

static void search_init(struct recv_fit_search *search, ssize_t nc, ssize_t dim,
			ssize_t ne)
{
	vector_init(&search->vector, dim * nc + ne);
}

static void search_add_constrs(struct recv_fit_search *search, ssize_t n)
{
	assert(n >= 0);

	ssize_t m0 = vector_dim(&search->vector);
	ssize_t m = m0 + n;
	vector_reinit(&search->vector, m);
}

static void search_deinit(struct recv_fit_search *search)
{
	vector_deinit(&search->vector);
}

static void search_set(struct recv_fit_search *search,
		       const struct recv_fit_resid *resid,
		       struct recv_fit_kkt *kkt)
{
	assert(!kkt->factored);

	/* determine the search direction */
	struct vector *s = &search->vector;
	vector_assign_copy(s, &resid->vector);
	vector_scale(s, -1.0);

	struct matrix smat = matrix_make(s, vector_dim(s), 1);
	ssize_t info = ldlfac_solve(&kkt->ldl, kkt->uplo,
				    &kkt->matrix, &smat);
	kkt->factored = true;

	assert(info == 0);
	(void)info;		// compile warning;
}

static void rgrad_init(struct recv_fit_rgrad *rgrad, ssize_t nc, ssize_t dim,
		       ssize_t ne)
{
	vector_init(&rgrad->vector, dim * nc + ne);
}

static void rgrad_add_constrs(struct recv_fit_rgrad *rgrad, ssize_t n)
{
	assert(n >= 0);

	ssize_t m0 = vector_dim(&rgrad->vector);
	ssize_t m = m0 + n;
	vector_reinit(&rgrad->vector, m);
}

static void rgrad_deinit(struct recv_fit_rgrad *rgrad)
{
	vector_deinit(&rgrad->vector);
}

static void rgrad_set(struct recv_fit_rgrad *rgrad,
		      const struct recv_fit_kkt *kkt,
		      const struct recv_fit_resid *resid)
{
	assert(!kkt->factored);

	struct vector *g = &rgrad->vector;
	matrix_mul(2.0, BLAS_NOTRANS, &kkt->matrix, &resid->vector, 0.0, g);
}

static enum recv_fit_task primal_dual_step(struct recv_fit *fit)
{
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
	search_set(&fit->search, &fit->prev->resid, &fit->kkt);
	ssize_t ce = recv_fit_constr_count(fit);
	ssize_t nparam = vector_dim(&fit->search.vector);
	struct vector coef_search =
	    vector_slice(&fit->search.vector, 0, nparam - ce);
	double smax = vector_max_abs(&coef_search);

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
	ssize_t it = 0;
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
		eval_step(fit->cur, &fit->constr,
			  &fit->prev->params, fit->step, &fit->search.vector,
			  fit->xmsgs, fit->ymsgs, fit->frame, &fit->model);

		if (!fit->cur->in_domain) {
			goto domain_error_f;
		}
		// compute the kkt matrix and residual gradient at the new point
		kkt_set(&fit->kkt, &fit->cur->loglik, &fit->constr);
		rgrad_set(&fit->rgrad, &fit->kkt, &fit->cur->resid);

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
		   size_t ncohort,
		   const size_t *cohorts, const struct recv_fit_ctrl *ctrl)
{
	assert(fit);
	assert(xmsgs);
	assert(f);
	assert(messages_max_to(xmsgs) < frame_send_count(f));
	assert(messages_max_from(xmsgs) < frame_recv_count(f));
	assert(!ymsgs || (messages_max_to(ymsgs) < frame_send_count(f)));
	assert(!ymsgs || (messages_max_from(ymsgs) < frame_recv_count(f)));
	assert(!ctrl || recv_fit_ctrl_valid(ctrl));

	const struct design *design = frame_recv_design(f);
	size_t dim = design_dim(design);
	size_t nc = ncohort;
	size_t ne = 0;

	if (!ctrl) {
		fit->ctrl = RECV_FIT_CTRL0;
	} else {
		fit->ctrl = *ctrl;
	}

	frame_clear(f);
	fit->frame = f;
	fit->xmsgs = xmsgs;
	fit->ymsgs = ymsgs ? ymsgs : xmsgs;

	recv_model_init(&fit->model, f, ncohort, cohorts, NULL);
	constr_init(&fit->constr);
	eval_init(&fit->eval[0], &fit->model, ne);
	eval_init(&fit->eval[1], &fit->model, ne);
	fit->prev = &fit->eval[0];
	fit->cur = &fit->eval[1];
	kkt_init(&fit->kkt, nc, dim, ne);
	search_init(&fit->search, nc, dim, ne);
	rgrad_init(&fit->rgrad, nc, dim, ne);

	// evaluate the model at zero
	eval_set(fit->cur, &fit->constr, NULL, NULL, fit->xmsgs, fit->ymsgs,
		 fit->frame, &fit->model);

	fit->dev0 = recv_fit_dev(fit);
}

#if 0
ssize_t recv_fit_add_constr_identify(struct recv_fit *fit)
{
	ssize_t i, dim = design_recv_dim(fit->design);
	ssize_t ic, nc = actors_cohort_count(fit->senders);
	ssize_t ne = recv_fit_constr_count(fit);
	const struct recv_loglik *ll = &fit->cur->loglik;
	// ssize_t nmsg = recv_loglik_count_sum(ll);
	//const struct matrix *ce = &fit->constr.ce;
	struct vector scale;
	struct matrix imatc;
	struct symeig eig;
	struct vector evals;
	struct matrix *nullspaces;
	struct matrix proj;
	ssize_t nulldim;

	vector_init(&scale, dim * nc);
	matrix_init(&imatc, dim, dim);
	symeig_init(&eig, dim, EIG_VEC);
	vector_init(&evals, dim);
	nullspaces = xcalloc(nc, sizeof(nullspaces[0]));

	nulldim = 0;

	for (ic = 0; ic < nc; ic++) {
		const struct recv_loglik_info *ll_info =
		    recv_loglik_info(ll, ic);
		// ssize_t nmsgc = recv_loglik_count(ll, ic);

		struct vector scalec = vector_slice(&scale, ic * dim, dim);
		matrix_assign_copy(&imatc, BLAS_NOTRANS, &ll_info->imat);

		// compute the scale
		for (i = 0; i < dim; i++) {
			double s2 = matrix_item(&imatc, i, i);
			double s;

			if (s2 < fit->ctrl.vartol) {
				s = 1.0;
			} else {
				s = sqrt(s2);
			}

			vector_set_item(&scalec, i, 1.0 / s);
		}

		// scale the information matrix
		matrix_scale_rows(&imatc, &scalec);
		matrix_scale_cols(&imatc, &scalec);

		// compute the eigendecomposition
		symeig_factor(&eig, BLAS_LOWER, &imatc, &evals);

		// compute the rank of nullspace
		ssize_t nulldimc = 0;
		for (nulldimc = 0; nulldimc < dim; nulldimc++) {
			if (vector_item(&evals, nulldimc) > fit->ctrl.eigtol)
				break;
		}

		// update the global nullspace dim
		nulldim += nulldimc;

		// copy the nullspace
		struct matrix nullc = matrix_slice_cols(&imatc, 0, nulldimc);
		matrix_init_copy(&nullspaces[ic], BLAS_NOTRANS, &nullc);

		// change back to the original scale
		matrix_scale_rows(&nullspaces[ic], &scalec);
	}

	// compute the projection of the constraints on to the nullsapce
	matrix_init(&proj, ne, nulldim);
	ssize_t off = 0;
	for (ic = 0; ic < nc; ic++) {
		struct matrix cec = matrix_slice_rows(ce, ic * dim, dim);

		struct matrix *nullc = &nullspaces[ic];
		ssize_t nulldimc = matrix_ncol(nullc);

		struct matrix projc = matrix_slice_cols(&proj, off, nulldimc);
		matrix_matmul(1.0, BLAS_TRANS, &cec, nullc, 0.0, &projc);

		off += nulldimc;
	}
	assert(off == nulldim);

	// compute the nullspace of the coefficient matrix
	struct matrix ncoefs;

	if (ne > 0) {
		struct vector s;
		struct matrix u, vt;
		struct svdfac svd;
		bool ok;

		vector_init(&s, MIN(ne, nulldim));
		matrix_init(&u, ne, ne);
		matrix_init(&vt, nulldim, nulldim);
		svdfac_init(&svd, ne, nulldim, SVD_ALL);
		ok = svdfac_factor(&svd, &proj, &s, &u, &vt);
		assert(ok);

		ssize_t rank, maxrank = MIN(ne, nulldim);
		for (rank = maxrank; rank > 0; rank--) {
			if (vector_item(&s, rank - 1) > fit->ctrl.svtol)
				break;
		}

		struct matrix vt0 =
		    matrix_slice_rows(&vt, rank, nulldim - rank);
		matrix_init_copy(&ncoefs, BLAS_TRANS, &vt0);

		svdfac_deinit(&svd);
		matrix_deinit(&vt);
		matrix_deinit(&u);
		vector_deinit(&s);

	} else {
		matrix_init(&ncoefs, nulldim, nulldim);
		matrix_assign_identity(&ncoefs);
	}

	// add the new constraints
	assert(matrix_nrow(&ncoefs) == nulldim);
	struct vector ce1;
	double be1 = 0.0;
	vector_init(&ce1, dim * nc);

	ssize_t ic1, nc1 = matrix_ncol(&ncoefs);
	for (ic1 = 0; ic1 < nc1; ic1++) {
		struct vector y = matrix_col(&ncoefs, ic1);

		vector_fill(&ce1, 0.0);
		off = 0;
		for (ic = 0; ic < nc; ic++) {
			struct matrix *nullc = &nullspaces[ic];
			ssize_t nulldimc = matrix_ncol(nullc);

			struct vector ce1c = vector_slice(&ce1, ic * dim, dim);
			struct vector yc = vector_slice(&y, off, nulldimc);
			matrix_mul(1.0, BLAS_NOTRANS, nullc, &yc, 0.0, &ce1c);
			off += nulldimc;
		}
		assert(off == nulldim);

		const char *fmt = "Ident%" SSIZE_FMT "";
		char zero;
		size_t len = snprintf(&zero, 1, fmt, ic + 1) + 1;
		char *name = xcalloc(len, sizeof(*name));
		snprintf(name, len, fmt, ic1 + 1);
		recv_fit_add_constr(fit, &ce1, be1, name);
		free(name);
	}

	vector_deinit(&ce1);
	matrix_deinit(&ncoefs);
	matrix_deinit(&proj);

	for (ic = 0; ic < nc; ic++) {
		matrix_deinit(&nullspaces[ic]);
	}
	free(nullspaces);

	vector_deinit(&evals);
	symeig_deinit(&eig);
	matrix_deinit(&imatc);
	vector_deinit(&scale);

	return nc1;
}
#endif

enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct matrix *coefs0)
{
	if (coefs0 != NULL) {
		eval_set(fit->cur, &fit->constr,
			 coefs0, NULL, fit->xmsgs, fit->ymsgs, fit->frame,
			 &fit->model);
	}
	kkt_set(&fit->kkt, &fit->cur->loglik, &fit->constr);
	fit->step = NAN;

	// determine initial value and gradient
	double f0 = fit->cur->resid.norm2;
	//double g0 = -2 * f0;

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
			 const double **weightsp, const size_t **indp,
			 size_t *nzp, double *valp, const char **namep)
{
	constr_get(&fit->constr, i, weightsp, indp, nzp, valp, namep);
}

static void _recv_fit_add_constrs(struct recv_fit *fit, ssize_t n)
{
	eval_add_constrs(&fit->eval[0], n);
	eval_add_constrs(&fit->eval[1], n);
	kkt_add_constrs(&fit->kkt, n);
	search_add_constrs(&fit->search, n);
	rgrad_add_constrs(&fit->rgrad, n);

}

void recv_fit_add_constr(struct recv_fit *fit, const double *weights,
			 const size_t *ind, size_t nz, double val,
			 const char *name)
{
	constr_add(&fit->constr, weights, ind, nz, val, name);
	_recv_fit_add_constrs(fit, 1);
}

void recv_fit_add_constr_set(struct recv_fit *fit, size_t i, size_t c,
			     double val)
{
	const struct design *d = frame_recv_design(fit->frame);
	size_t dim = design_dim(d);
	size_t nc = recv_model_cohort_count(&fit->model);

	constr_add_set(&fit->constr, dim, nc, i, c, val);
	_recv_fit_add_constrs(fit, 1);
}

void recv_fit_add_constr_eq(struct recv_fit *fit, size_t i1, size_t c1,
			    size_t i2, size_t c2)
{
	const struct design *d = frame_recv_design(fit->frame);
	size_t dim = design_dim(d);
	size_t nc = recv_model_cohort_count(&fit->model);

	constr_add_eq(&fit->constr, dim, nc, i1, c1, i2, c2);
	_recv_fit_add_constrs(fit, 1);
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
	double dev = 0.0;
	ssize_t ic, nc = recv_model_cohort_count(&fit->model);
	for (ic = 0; ic < nc; ic++) {
		const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
		ssize_t nrecv = info->nrecv;
		double avg_dev = info->dev;
		dev += nrecv * avg_dev;
	}

	return dev;
}

const struct matrix *recv_fit_coefs(const struct recv_fit *fit)
{
	assert(fit);
	return &fit->cur->coefs;
}

const struct vector *recv_fit_duals(const struct recv_fit *fit)
{
	assert(fit);
	return &fit->cur->duals;
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
