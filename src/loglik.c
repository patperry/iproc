#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "sloglik.h"
#include "loglik.h"


void recv_loglik_init(struct recv_loglik *ll,
		      struct model *m)
{
	assert(ll);
	assert(m);

	const struct design *design = model_design(m);
	ssize_t isend, nsend = design_send_count(design);
	
	ll->model = model_ref(m);
	
	array_init(&ll->slogliks, sizeof(struct recv_sloglik));
	array_set_capacity(&ll->slogliks, nsend);
	for (isend = 0; isend < nsend; isend++) {
		struct recv_sloglik *sll = array_add(&ll->slogliks, NULL);
		recv_sloglik_init(sll, m, isend);
	}

	vector_init(&ll->grad, design_recv_dim(design));
	ll->grad_cached = true;
	ll->nsend = 0;
	ll->nrecv = 0;
	ll->last = NULL;
}

void recv_loglik_deinit(struct recv_loglik *ll)
{
	assert(ll);

	struct recv_sloglik *sll;
	ARRAY_FOREACH(sll, &ll->slogliks) {
		recv_sloglik_deinit(sll);
	}
	array_deinit(&ll->slogliks);
	vector_deinit(&ll->grad);
	model_free(ll->model);
}

struct recv_loglik *recv_loglik_alloc(struct model *m, struct messages *msgs)
{
	struct recv_loglik *ll = xcalloc(1, sizeof(*ll));
	
	struct frame f;
	
	frame_init(&f, model_design(m));
	recv_loglik_init(ll, m);
	recv_loglik_add_all(ll, &f, msgs);
	frame_deinit(&f);
	
	return ll;
}

void recv_loglik_free(struct recv_loglik *ll)
{
	if (ll) {
		recv_loglik_deinit(ll);
		xfree(ll);
	}
}

void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f,
		     const struct message *msg)
{
	struct recv_sloglik *sll = array_item(&ll->slogliks, msg->from);
	recv_sloglik_add(sll, f, msg->to, msg->nto);
	ll->grad_cached = false;
	ll->nsend += 1;
	ll->nrecv += msg->nto;
	ll->last = sll;
}

void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f,
			 const struct messages *msgs)
{
	struct messages_iter it;
	const struct message *msg;
	ssize_t i, n;
	
	struct model *m = ll->model;

	MESSAGES_FOREACH(it, msgs) {
		double t = MESSAGES_TIME(it);
		
		while (frame_next_change(f) <= t) {
			model_update(m, f);
			frame_advance(f);			
		}
		if (frame_time(f) < t) {
			model_update(m, f);
			frame_advance_to(f, t);			
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			recv_loglik_add(ll, f, msg);
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(f, msg);
		}
	}
	
}
			 

double recv_loglik_value(const struct recv_loglik *ll)
{
	assert(ll);
	
	double value = 0.0;
	struct recv_sloglik *sll;

	ARRAY_FOREACH(sll, &ll->slogliks) {
		value += recv_sloglik_value(sll);
	}

	return value;
}

static void recv_loglik_axpy_grad_nocache(double alpha,
					  const struct recv_loglik *ll,
					  struct vector *y)
{
	assert(ll);
	assert(y);
	assert(vector_dim(y) == model_dim(ll->model));
	
	struct recv_sloglik *sll;
	
	ARRAY_FOREACH(sll, &ll->slogliks) {
		recv_sloglik_axpy_grad(alpha, sll, y);
	}
}

static void recv_loglik_cache_grad(struct recv_loglik *ll)
{
	assert(ll);

	vector_fill(&ll->grad, 0.0);
	recv_loglik_axpy_grad_nocache(1.0, ll, &ll->grad);
	ll->grad_cached = true;
}

void recv_loglik_axpy_grad(double alpha,
			   const struct recv_loglik * ll,
			   struct vector *y)
{
	assert(ll);
	assert(y);
	assert(vector_dim(y) == model_dim(ll->model));
	
	if (!ll->grad_cached) {
		recv_loglik_cache_grad((struct recv_loglik *)ll);
	}
	
	vector_axpy(alpha, &ll->grad, y);
}

ssize_t recv_loglik_count(const struct recv_loglik *ll)
{
	assert(ll);
	return ll->nrecv;
}

double recv_loglik_avg_dev(const struct recv_loglik *ll)
{
	assert(ll);
	double dev, dev_avg;
	ssize_t ntot, n;
	struct recv_sloglik *sll;
	
	ntot = 0;
	dev_avg = 0.0;
	
	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			dev = recv_sloglik_avg_dev(sll);
			dev_avg += n * (dev - dev_avg) / (n + ntot);
			ntot += n;
		}
	}
	assert(ntot == recv_loglik_count(ll));
	
	return dev_avg;
}

double recv_loglik_last_dev(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);
	
	return recv_sloglik_last_dev(ll->last);
}

void recv_loglik_axpy_avg_mean(double alpha, const struct recv_loglik *ll, struct vector *y)
{
	struct vector avg_mean, diff;
	vector_init(&avg_mean, model_dim(ll->model));	
	vector_init(&diff, model_dim(ll->model));
	ssize_t ntot, n;
	struct recv_sloglik *sll;
	
	ntot = 0;
	
	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			ntot += n;			
			vector_assign_copy(&diff, &avg_mean);
			recv_sloglik_axpy_avg_mean(-1.0, sll, &diff);
			vector_axpy(-((double)n)/ntot, &diff, &avg_mean);
		}
	}
	assert(ntot == recv_loglik_count(ll));
	
	vector_axpy(alpha, &avg_mean, y);
	vector_deinit(&diff);
	vector_deinit(&avg_mean);
}

void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll, struct vector *y)
{
	assert(ll);
	assert(ll->last);
	recv_sloglik_axpy_last_mean(alpha, ll->last, y);
}

void recv_loglik_axpy_avg_imat(double alpha, const struct recv_loglik *ll, struct matrix *y)
{
	struct matrix avg_imat, diff;
	ssize_t dim = model_dim(ll->model);
	matrix_init(&avg_imat, dim, dim);
	matrix_init(&diff, dim, dim);
	ssize_t ntot, n;
	struct recv_sloglik *sll;
	
	ntot = 0;
	
	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			ntot += n;			
			matrix_assign_copy(&diff, &avg_imat);
			recv_sloglik_axpy_avg_imat(-1.0, sll, &diff);
			matrix_axpy(-((double)n)/ntot, &diff, &avg_imat);
		}
	}
	assert(ntot == recv_loglik_count(ll));
	
	matrix_axpy(alpha, &avg_imat, y);
	matrix_deinit(&diff);
	matrix_deinit(&avg_imat);
}

void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll, struct matrix *y)
{
	assert(ll);
	assert(ll->last);
	recv_sloglik_axpy_last_imat(alpha, ll->last, y);
}

