#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linesearch.h"
#include "fit.h"

static void
eval_objective(struct recv_loglik *loglik,
	       double penalty, double *valuep, struct vector *grad)
{
	struct model *model = loglik->model;
	struct vector *coefs = &model->coefs;

	double n = loglik->nrecv;
	double ll_value = -recv_loglik_avg_dev(loglik) / 2.0;

	double norm = vector_norm(coefs);
	double norm2 = norm * norm;
	double value = -ll_value / n + 0.5 * (penalty / n) * norm2;
	*valuep = value;

	vector_assign_copy(grad, coefs);
	vector_scale(grad, penalty / n);
	recv_loglik_axpy_avg_score(1.0, loglik, grad);
}

static void iproc_fit_init(iproc_fit * fit)
{
	ssize_t dim = model_dim(fit->model);

	struct bfgs_ctrl ctrl = BFGS_CTRL0;
	bfgs_init(&fit->opt, dim, &ctrl);
	vector_init(&fit->grad, dim);

	fit->loglik = recv_loglik_alloc(fit->model, fit->messages);

	eval_objective(fit->loglik, fit->penalty, &fit->f, &fit->grad);

	fit->xnext =
	    bfgs_start(&fit->opt, &fit->model->coefs, fit->f, &fit->grad);
}

iproc_fit *iproc_fit_new(struct model *model0,
			 struct messages *messages, double penalty)
{
	assert(model0);
	assert(messages);
	assert(penalty >= 0.0);
	assert(messages_max_from(messages) < model_sender_count(model0));
	assert(messages_max_to(messages) < model_receiver_count(model0));

	iproc_fit *fit = xcalloc(1, sizeof(*fit));

	fit->model = model_ref(model0);
	fit->messages = messages_ref(messages);
	fit->penalty = penalty;
	fit->loglik = NULL;

	iproc_fit_init(fit);
	return fit;
}

void iproc_fit_free(iproc_fit * fit)
{
	if (fit) {
		recv_loglik_free(fit->loglik);
		bfgs_deinit(&fit->opt);
		vector_deinit(&fit->grad);
		messages_free(fit->messages);
		model_free(fit->model);
		xfree(fit);
	}
}

void iproc_fit_step(iproc_fit * fit)
{
	assert(fit);
	assert(fit->xnext);

	struct model *model = fit->model;
	struct messages *messages = fit->messages;
	double penalty = fit->penalty;
	struct design *design = design_ref(model->design);
	struct recv_loglik *loglik = fit->loglik;

	model_free(model);
	model = model_alloc(design, fit->xnext);

	/* Update the loglik, value, and gradient */
	recv_loglik_free(loglik);
	loglik = recv_loglik_alloc(model, messages);
	eval_objective(loglik, penalty, &fit->f, &fit->grad);

	fit->xnext = bfgs_advance(&fit->opt, fit->f, &fit->grad);
}

bool iproc_fit_converged(iproc_fit * fit)
{
	return bfgs_converged(&fit->opt);
}
