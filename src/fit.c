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
	assert(coefs);
	frame_clear(&fit->frame);
	model_clear(&fit->model);
	
	// update coefs
	model_set(&fit->model, &fit->frame, coefs);
	coefs = model_recv_coefs(&fit->model);
	
	recv_loglik_clear(&fit->loglik);
	recv_loglik_add_all(&fit->loglik, &fit->frame, fit->msgs);
	
	double n = recv_loglik_count(&fit->loglik);
	double dev = recv_loglik_avg_dev(&fit->loglik);
	double norm = vector_norm(coefs);
	double norm2 = norm * norm;
	
	fit->f = 0.5 * (dev + (fit->penalty / n) * norm2);

	vector_assign_copy(&fit->grad, coefs);
	vector_scale(&fit->grad, fit->penalty / n);
	recv_loglik_axpy_avg_score(-1.0, &fit->loglik, &fit->grad);

	assert(isfinite(fit->f));
	assert(isfinite(vector_norm(&fit->grad)));
	
	/* hessian */
	matrix_fill(&fit->imat_evec, 0.0);
	recv_loglik_axpy_avg_imat(1.0, &fit->loglik, &fit->imat_evec);
	ssize_t i, dim = vector_dim(coefs);
	for (i = 0; i < dim; i++) {
		double *x = matrix_item_ptr(&fit->imat_evec, i, i);
		*x += (fit->penalty / n);
	}
	
	enum matrix_uplo uplo = UPLO_LOWER; // UPLO_UPPER would also be ok
	bool ok = symeig_factor(&fit->eig, uplo, &fit->imat_evec, &fit->imat_eval);
	double emax = dim > 0 ? vector_item(&fit->imat_eval, dim - 1) : 0;
	assert(ok);
	
	for (i = 0; i < dim; i++) {
		double eval = vector_item(&fit->imat_eval, i);
		assert(eval > -1e-5 * MAX(1, emax));
		if (eval >= 1e-5)
			break;
	}
	fit->imat_rank = dim - i;
	assert(fit->imat_rank > 0);
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
	
	vector_init(&fit->grad, dim);
	matrix_init(&fit->imat_evec, dim, dim);
	vector_init(&fit->imat_eval, dim);
	symeig_init(&fit->eig, dim, EIG_VEC);

	if (!coefs0) {
		vector_fill(&fit->grad, 0.0); // use grad for temp storage
		coefs0 = &fit->grad;
	}
	
	recv_fit_set(fit, coefs0);
	
	struct linesearch_ctrl lsctrl = LINESEARCH_CTRL0;
	fit->lsctrl = lsctrl;
}

void recv_fit_deinit(struct recv_fit *fit)
{
	assert(fit);
	symeig_deinit(&fit->eig);
	vector_deinit(&fit->imat_eval);
	matrix_deinit(&fit->imat_evec);
	vector_deinit(&fit->grad);
	recv_loglik_deinit(&fit->loglik);
	model_deinit(&fit->model);
	frame_deinit(&fit->frame);
}

bool recv_fit_step(struct recv_fit *fit)
{
	assert(fit);

	ssize_t dim = design_recv_dim(fit->design);
	ssize_t rank = fit->imat_rank;
	ssize_t off = dim - rank;
	struct matrix U = matrix_slice_cols(&fit->imat_evec, off, rank);
	struct vector l = vector_slice(&fit->imat_eval, off, rank);
	struct vector Ut_g, search;
	vector_init(&Ut_g, rank);
	vector_init(&search, dim);
	
	matrix_mul(-1.0, TRANS_TRANS, &U, &fit->grad, 0.0, &Ut_g);
	vector_div(&Ut_g, &l);
	matrix_mul(+1.0, TRANS_NOTRANS, &U, &Ut_g, 0.0, &search);
	
	struct vector x0;
	vector_init_copy(&x0, model_recv_coefs(&fit->model));
	double f0 = fit->f;
	double g0 = vector_dot(&search, &fit->grad);
	
	assert(g0 < 0);
	
	printf("[[[[ g0 = %.22f ]]]]\n", g0);
	
	if (-g0 <= 0)
		return false; // converged
	
	struct vector x;
	vector_init(&x, dim);
		
	double stp, f, g;
	
	stp = MIN(1.0, 1.0 / vector_max_abs(&search));
	
	
	linesearch_start(&fit->ls, stp, f0, g0, &fit->lsctrl);
	enum linesearch_task task;
		
	do {
		stp = linesearch_step(&fit->ls);
		vector_assign_copy(&x, &x0);
		vector_axpy(stp, &search, &x);
		recv_fit_set(fit, &x);
		f = fit->f;
		g = vector_dot(&search, &fit->grad);
		task = linesearch_advance(&fit->ls, f, g);
	} while (task == LINESEARCH_STEP && !linesearch_sdec(&fit->ls));
	
	printf("[[[[ stp = %.22f ]]]]\n", stp);
	
	return (task == LINESEARCH_CONV);
}

bool recv_fit_converged(const struct recv_fit *fit)
{
	return false;
}

const char *recv_fit_errmsg(const struct recv_fit *fit)
{
	return NULL;
}
