#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "linesearch.h"
#include "fit.h"


static void recv_fit_set(struct recv_fit *fit, const struct vector *coefs)
{
	frame_clear(&fit->frame);
	model_clear(&fit->model);
	
	// update coefs
	coefs = model_recv_coefs(&fit->model);
	
	recv_loglik_clear(&fit->loglik);
	recv_loglik_add_all(&fit->loglik, &fit->frame, fit->msgs);
	
	double n = recv_loglik_count(&fit->loglik);
	double ll = -0.5 * recv_loglik_avg_dev(&fit->loglik);
	double norm = vector_norm(coefs);
	double norm2 = norm * norm;
	
	fit->f = ll + 0.5 * fit->penalty * (norm2 / n);

	vector_assign_copy(&fit->grad, coefs);
	vector_scale(&fit->grad, fit->penalty / n);
	recv_loglik_axpy_avg_score(-1.0, &fit->loglik, &fit->grad);
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
	model_init(&fit->model, fit->design, coefs0);
	recv_loglik_init(&fit->loglik, &fit->model);
	
	struct bfgs_ctrl ctrl = BFGS_CTRL0;
	bfgs_init(&fit->opt, dim, &ctrl);
	vector_init(&fit->grad, dim);

	recv_fit_set(fit, coefs0);
	fit->xnext = bfgs_start(&fit->opt, model_recv_coefs(&fit->model), fit->f, &fit->grad);
}

void recv_fit_deinit(struct recv_fit *fit)
{
	assert(fit);
	vector_deinit(&fit->grad);
	bfgs_deinit(&fit->opt);	
	recv_loglik_deinit(&fit->loglik);
	model_deinit(&fit->model);
	frame_deinit(&fit->frame);
}

void recv_fit_step(struct recv_fit *fit)
{
	assert(fit);
	assert(fit->xnext);

	recv_fit_set(fit, fit->xnext);
	fit->xnext = bfgs_advance(&fit->opt, fit->f, &fit->grad);
}

bool recv_fit_converged(struct recv_fit *fit)
{
	return bfgs_converged(&fit->opt);
}
