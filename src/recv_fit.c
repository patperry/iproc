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




static void _eval_update(struct recv_fit_eval *eval,
			 const struct messages *xmsgs,
			 const struct messages *ymsgs,
			 struct frame *frame, struct recv_model *model,
			 int moments)
{
	// clear frame; set model coefs
	struct history *h = frame_history(frame);
	history_clear(h);
	recv_model_set_par(model, &eval->params.params);
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
			recv_loglik_add(&eval->loglik, msg);
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

static void set_dev0(struct recv_fit *fit)
{
	const struct frame *f = fit->frame;
	size_t nrecv = frame_recv_count(f);

	if (!frame_has_loops(f)) {
		assert(nrecv > 0);
		nrecv--;
	}

	fit->dev0 = 2 * messages_recv_count(fit->ymsgs) * log(nrecv);

	//eval_set(fit->cur, NULL, fit->xmsgs, fit->ymsgs, fit->frame,
	//	   &fit->model, 2);
	//fit->dev0 = recv_fit_dev(fit);
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

	set_dev0(fit);
}


enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct recv_fit_params *params0)
{
	eval_set(fit->cur, params0, fit->xmsgs, fit->ymsgs, fit->frame,
		 &fit->model, 2);
	resid_set(&fit->cur->resid, &fit->cur->loglik, fit->constr,
		  &fit->cur->params);
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


size_t constr_add_identify_recv_fit(struct constr *c, struct frame *f,
				    const struct messages *xmsgs,
				    const struct messages *ymsgs)
{
	struct recv_model model;
	struct recv_fit_eval eval;

	recv_model_init(&model, f, NULL);
	eval_init(&eval, &model, c);
	eval_set(&eval, NULL, xmsgs, ymsgs, f, &model, 2);

	size_t n = recv_coefs_dim(f);
	double *imatp = xcalloc(n * (n + 1) / 2, sizeof(double));

	recv_loglik_axpy_imat(1.0, &eval.loglik, imatp);

	size_t nadd = constr_add_identify(c, imatp, RECV_LOGLIK_UPLO);

	free(imatp);
	eval_deinit(&eval);
	recv_model_deinit(&model);

	return nadd;
}


enum recv_fit_task recv_fit_advance(struct recv_fit *fit)
{
	assert(fit);
	assert(fit->task == RECV_FIT_STEP);

	fit->task = primal_dual_step(fit);
	return fit->task;
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
