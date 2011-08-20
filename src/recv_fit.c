#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linesearch.h"
#include "util.h"
#include "recv_fit.h"

static void constr_init(struct recv_fit_constr *constr, ssize_t nc, ssize_t dim)
{
	constr->dim = dim;
	constr->nc = nc;
	matrix_init(&constr->ce, nc * dim, 0);
	vector_init(&constr->be, 0);
}

static void constr_deinit(struct recv_fit_constr *constr)
{
	vector_deinit(&constr->be);
	matrix_deinit(&constr->ce);
}

static ssize_t constr_count(const struct recv_fit_constr *constr)
{
	assert(matrix_ncol(&constr->ce) == vector_dim(&constr->be));
	return vector_dim(&constr->be);
}

static void constr_add(struct recv_fit_constr *constr, const struct vector *ce,
		       double be)
{
	ssize_t ne0 = vector_dim(&constr->be);
	ssize_t ne = ne0 + 1;

	matrix_reinit(&constr->ce, matrix_nrow(&constr->ce), ne);

	struct vector ce1 = matrix_col(&constr->ce, ne0);
	vector_assign_copy(&ce1, ce);

	vector_reinit(&constr->be, ne);
	vector_set_item(&constr->be, ne0, be);
}

static void constr_add_set(struct recv_fit_constr *constr, ssize_t i, ssize_t c,
			   double val)
{
	assert(constr);
	assert(0 <= i && i < constr->dim);
	assert(0 <= c && c < constr->nc);
	assert(isfinite(val));

	ssize_t ne0 = vector_dim(&constr->be);
	ssize_t ne = ne0 + 1;

	matrix_reinit(&constr->ce, matrix_nrow(&constr->ce), ne);
	matrix_set_item(&constr->ce, i + c * constr->dim, ne0, 1.0);

	vector_reinit(&constr->be, ne);
	vector_set_item(&constr->be, ne0, val);
}

static void constr_add_eq(struct recv_fit_constr *constr, ssize_t i1,
			  ssize_t c1, ssize_t i2, ssize_t c2)
{
	assert(constr);
	assert(0 <= i1 && i1 < constr->dim);
	assert(0 <= i2 && i2 < constr->dim);
	assert(0 <= c1 && c1 < constr->nc);
	assert(0 <= c2 && c2 < constr->nc);
	assert(!(i1 == i2 && c1 == c2));

	ssize_t ne0 = vector_dim(&constr->be);
	ssize_t ne = ne0 + 1;

	matrix_reinit(&constr->ce, matrix_nrow(&constr->ce), ne);
	matrix_set_item(&constr->ce, i1 + c1 * constr->dim, ne0, +1.0);
	matrix_set_item(&constr->ce, i2 + c2 * constr->dim, ne0, -1.0);

	vector_reinit(&constr->be, ne);
	vector_set_item(&constr->be, ne0, 0.0);
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
		      const struct recv_fit_constr *constr,
		      const struct vector *params)
{
	struct vector *r = &resid->vector;
	const struct recv_model *m = recv_loglik_model(ll);
	ssize_t ntot = recv_loglik_count_sum(ll);
	ssize_t dim = recv_model_dim(m);
	ssize_t ne = constr_count(constr);
	ssize_t ic, nc = recv_model_cohort_count(m);

	/* r1 is the dual residual: grad(f) + ce * nu,
	 *   where grad(f) = grad(nll)
	 *                 = -(score)
	 */
	struct vector r1 = vector_slice(r, 0, nc * dim);
	struct vector duals = vector_slice(params, nc * dim, ne);
	matrix_mul(1.0, TRANS_NOTRANS, &constr->ce, &duals, 0.0, &r1);

	for (ic = 0; ic < nc; ic++) {
		const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
		const struct vector *score = &info->score;
		ssize_t n = recv_loglik_count(ll, ic);

		struct vector r1c = vector_slice(&r1, ic * dim, dim);
		vector_axpy(-((double)n) / ntot, score, &r1c);
	}

	// r2 is the primal residual: ce' * x - be
	struct vector r2 = vector_slice(r, nc * dim, ne);
	struct vector primals = vector_slice(params, 0, nc * dim);

	vector_assign_copy(&r2, &constr->be);
	matrix_mul(1.0, TRANS_TRANS, &constr->ce, &primals, -1.0, &r2);

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
			 const struct recv_fit_constr *constr,
			 const struct messages *msgs, struct frame *frame,
			 struct recv_model *model)
{
	// clear frame; set model coefs
	frame_clear(frame);
	recv_model_set_coefs(model, &eval->coefs);

	// set loglik
	recv_loglik_clear(&eval->loglik);
	recv_loglik_add_all(&eval->loglik, frame, msgs);

	// set resid
	resid_set(&eval->resid, &eval->loglik, constr, &eval->params);

	// set in_domain
	eval->in_domain = isfinite(eval->resid.norm2);

}

static void eval_set(struct recv_fit_eval *eval,
		     const struct recv_fit_constr *constr,
		     const struct matrix *coefs, const struct vector *duals,
		     const struct messages *msgs, struct frame *frame,
		     struct recv_model *model)
{
	ssize_t dim = recv_model_dim(model);
	ssize_t nc = recv_model_cohort_count(model);
	ssize_t ne = constr_count(constr);

	struct vector pcoefs_data = vector_slice(&eval->params, 0, dim * nc);
	struct matrix pcoefs = matrix_make(&pcoefs_data, dim, nc);

	if (coefs) {
		matrix_assign_copy(&pcoefs, TRANS_NOTRANS, coefs);
	} else {
		matrix_fill(&pcoefs, 0.0);
	}

	struct vector pduals = vector_slice(&eval->params, dim * nc, ne);

	if (duals) {
		vector_assign_copy(&pduals, duals);
	} else {
		vector_fill(&pduals, 0.0);
	}

	_eval_update(eval, constr, msgs, frame, model);
}

static void eval_step(struct recv_fit_eval *eval,
		      const struct recv_fit_constr *constr,
		      const struct vector *params0, double scale,
		      const struct vector *dir, const struct messages *msgs,
		      struct frame *frame, struct recv_model *model)
{
	vector_assign_copy(&eval->params, params0);
	vector_axpy(scale, dir, &eval->params);
	_eval_update(eval, constr, msgs, frame, model);
}

static void kkt_init(struct recv_fit_kkt *kkt, ssize_t nc, ssize_t dim,
		     ssize_t ne)
{
	ssize_t n = nc * dim + ne;

	kkt->factored = false;
	kkt->uplo = UPLO_LOWER;
	matrix_init(&kkt->matrix, n, n);
	ldlfac_init(&kkt->ldl, n);
}

static void kkt_add_constrs(struct recv_fit_kkt *kkt, ssize_t n)
{
	assert(n >= 0);
	ssize_t m = matrix_nrow(&kkt->matrix) + n;

	kkt->factored = false;
	kkt->uplo = UPLO_LOWER;
	matrix_reinit(&kkt->matrix, m, m);
	ldlfac_reinit(&kkt->ldl, m);
}

static void kkt_deinit(struct recv_fit_kkt *kkt)
{
	ldlfac_deinit(&kkt->ldl);
	matrix_deinit(&kkt->matrix);
}

static void kkt_set(struct recv_fit_kkt *kkt, const struct recv_loglik *ll,
		    const struct recv_fit_constr *constr)
{
	struct matrix *k = &kkt->matrix;
	const struct recv_model *m = recv_loglik_model(ll);
	const struct actors *s = recv_model_senders(m);
	ssize_t dim = recv_model_dim(m);
	ssize_t ne = constr_count(constr);
	ssize_t ic, nc = actors_cohort_count(s);
	ssize_t ntot = recv_loglik_count_sum(ll);

	matrix_fill(k, 0.0);

	// k11
	for (ic = 0; ic < nc; ic++) {
		struct recv_loglik_info *info = recv_loglik_info(ll, ic);
		ssize_t n = recv_loglik_count(ll, ic);
		struct matrix h = matrix_slice(k, ic * dim, ic * dim, dim, dim);
		matrix_assign_copy(&h, TRANS_NOTRANS, &info->imat);
		matrix_scale(&h, ((double)n) / ntot);
	}

	struct matrix k21 = matrix_slice(k, nc * dim, 0, ne, nc * dim);
	matrix_assign_copy(&k21, TRANS_TRANS, &constr->ce);

	struct matrix k12 = matrix_slice(k, 0, nc * dim, nc * dim, ne);
	matrix_assign_copy(&k12, TRANS_NOTRANS, &constr->ce);

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
	(void)info; // compile warning;
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
	matrix_mul(2.0, TRANS_NOTRANS, &kkt->matrix, &resid->vector, 0.0, g);
}

static enum recv_fit_task primal_dual_step(struct recv_fit *fit)
{
	// swap prev and cur evals
	struct recv_fit_eval *tmp = fit->prev;
	fit->prev = fit->cur;
	fit->cur = tmp;

	// determine initial value and gradient
	double f0 = fit->prev->resid.norm2;
	//double g0 = -2 * f0;

	// stop early if residual is below tolerance
	if (f0 <= fit->ctrl.gtol * fit->ctrl.gtol)
		return RECV_FIT_CONV;

	// determine the search direction
	search_set(&fit->search, &fit->prev->resid, &fit->kkt);
	ssize_t ce = recv_fit_constr_count(fit);
	ssize_t nparam = vector_dim(&fit->search.vector);
	struct vector coef_search = vector_slice(&fit->search.vector, 0, nparam - ce);
	double smax = vector_max_abs(&coef_search);

	// set up the linesearch control parameters
	struct linesearch_ctrl ctrl = fit->ctrl.ls;
	ctrl.stpmax = MIN(ctrl.stpmax, 100 / MAX(1.0, smax));
	ctrl.stpmin = MIN(ctrl.stpmin, 1e-12 * ctrl.stpmax);
	double stp0 = MIN(1.0, 4.0 / MAX(1.0, smax));
	double stp = stp0;

	fprintf(stderr, "> smax = %.8f   stpmax = %.8f   stp0 = %.8f\n", smax, ctrl.stpmax, stp0);
	
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
			  fit->msgs, &fit->frame, &fit->model);

		if (!fit->cur->in_domain) {
			goto domain_error_f;
		}
		// compute the kkt matrix and residual gradient at the new point
		kkt_set(&fit->kkt, &fit->cur->loglik, &fit->constr);
		rgrad_set(&fit->rgrad, &fit->kkt, &fit->cur->resid);

		double f = fit->cur->resid.norm2;
		//double g = vector_dot(&fit->search.vector, &fit->rgrad.vector);

		//if (!isfinite(g)) {
		//	goto domain_error_g;
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
	} while (it < fit->ctrl.ls_maxit
		 && task == LINESEARCH_STEP); // && !linesearch_sdec(&fit->ls));

	if (task == LINESEARCH_WARN_XTOL) {
		return RECV_FIT_ERR_XTOL;
	} else if (task == LINESEARCH_CONV) {
		return RECV_FIT_STEP;
	} else {
		return RECV_FIT_ERR_LNSRCH;
	}
}

void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct actors *senders,
		   const struct recv_fit_ctrl *ctrl)
{
	assert(fit);
	assert(msgs);
	assert(design);
	assert(messages_max_to(msgs) < design_send_count(design));
	assert(messages_max_from(msgs) < design_recv_count(design));
	assert(!ctrl || recv_fit_ctrl_valid(ctrl));

	ssize_t dim = design_recv_dim(design);
	ssize_t nc = actors_cohort_count(senders);
	ssize_t ne = 0;

	if (!ctrl) {
		fit->ctrl = RECV_FIT_CTRL0;
	} else {
		fit->ctrl = *ctrl;
	}

	fit->msgs = msgs;
	fit->design = design;
	fit->senders = senders;

	frame_init(&fit->frame, fit->design);
	recv_model_init(&fit->model, &fit->frame, fit->senders, NULL);
	constr_init(&fit->constr, nc, dim);
	eval_init(&fit->eval[0], &fit->model, ne);
	eval_init(&fit->eval[1], &fit->model, ne);
	fit->prev = &fit->eval[0];
	fit->cur = &fit->eval[1];
	kkt_init(&fit->kkt, nc, dim, ne);
	search_init(&fit->search, nc, dim, ne);
	rgrad_init(&fit->rgrad, nc, dim, ne);
	
	// evaluate the model at zero
	eval_set(fit->cur, &fit->constr, NULL, NULL, fit->msgs, &fit->frame,
		 &fit->model);
}


ssize_t recv_fit_add_constr_identify(struct recv_fit *fit)
{
	ssize_t i, dim = design_recv_dim(fit->design);
	ssize_t ic, nc = actors_cohort_count(fit->senders);
	ssize_t ne = constr_count(&fit->constr);
	const struct recv_loglik *ll = &fit->cur->loglik;	
	// ssize_t nmsg = recv_loglik_count_sum(ll);
	const struct matrix *ce = &fit->constr.ce;
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
		const struct recv_loglik_info *ll_info = recv_loglik_info(ll, ic);
		// ssize_t nmsgc = recv_loglik_count(ll, ic);
		
		struct vector scalec = vector_slice(&scale, ic * dim, dim);
		matrix_assign_copy(&imatc, TRANS_NOTRANS, &ll_info->imat);
		
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
		symeig_factor(&eig, UPLO_LOWER, &imatc, &evals);
		
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
		matrix_init_copy(&nullspaces[ic], TRANS_NOTRANS, &nullc);
		
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
		matrix_matmul(1.0, TRANS_TRANS, &cec, nullc, 0.0, &projc);
		
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
		
		struct matrix vt0 = matrix_slice_rows(&vt, rank, nulldim - rank);
		matrix_init_copy(&ncoefs, TRANS_TRANS, &vt0);
		
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
			matrix_mul(1.0, TRANS_NOTRANS, nullc, &yc, 0.0, &ce1c);
			off += nulldimc;
		}
		assert(off == nulldim);
		
		recv_fit_add_constr(fit, &ce1, be1);
	}
	
	vector_deinit(&ce1);
	matrix_deinit(&ncoefs);
	matrix_deinit(&proj);

	for (ic = 0; ic < nc; ic++) {
		matrix_deinit(&nullspaces[ic]);
	}
	xfree(nullspaces);
	
	vector_deinit(&evals);
	symeig_deinit(&eig);
	matrix_deinit(&imatc);
	vector_deinit(&scale);
	
	return nc1;
}


enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct matrix *coefs0)
{
	if (coefs0 != NULL) {
		eval_set(fit->cur, &fit->constr, coefs0, NULL, fit->msgs, &fit->frame,
			 &fit->model);
	}
	kkt_set(&fit->kkt, &fit->cur->loglik, &fit->constr);
	fit->step = NAN;

	// determine initial value and gradient
	double f0 = fit->cur->resid.norm2;
	double g0 = -2 * f0;

	// stop early if residual is below tolerance
	if (-g0 < fit->ctrl.gtol * fit->ctrl.gtol) {
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
	frame_deinit(&fit->frame);
}

ssize_t recv_fit_constr_count(const struct recv_fit *fit)
{
	return constr_count(&fit->constr);
}

void recv_fit_get_constr(const struct recv_fit *fit, const struct matrix **ce,
			 const struct vector **be)
{
	*ce = &fit->constr.ce;
	*be = &fit->constr.be;
}

static void _recv_fit_add_constrs(struct recv_fit *fit, ssize_t n)
{
	eval_add_constrs(&fit->eval[0], n);
	eval_add_constrs(&fit->eval[1], n);
	kkt_add_constrs(&fit->kkt, n);
	search_add_constrs(&fit->search, n);
	rgrad_add_constrs(&fit->rgrad, n);

}

void recv_fit_add_constr(struct recv_fit *fit, const struct vector *ce,
			 double be)
{
	constr_add(&fit->constr, ce, be);
	_recv_fit_add_constrs(fit, 1);
}

void recv_fit_add_constr_set(struct recv_fit *fit, ssize_t i, ssize_t c,
			     double val)
{
	constr_add_set(&fit->constr, i, c, val);
	_recv_fit_add_constrs(fit, 1);
}

void recv_fit_add_constr_eq(struct recv_fit *fit, ssize_t i1, ssize_t c1,
			    ssize_t i2, ssize_t c2)
{
	constr_add_eq(&fit->constr, i1, c1, i2, c2);
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
