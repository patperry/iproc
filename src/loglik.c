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

double recv_loglik_mean_dev(const struct recv_loglik *ll)
{
	assert(ll);
	double dev, dev_mean;
	ssize_t ntot, n;
	struct recv_sloglik *sll;
	
	ntot = 0;
	dev_mean = 0.0;
	
	ARRAY_FOREACH(sll, &ll->slogliks) {
		n = recv_sloglik_count(sll);
		if (n > 0) {
			dev = recv_sloglik_mean_dev(sll);
			dev_mean += n * (dev - dev_mean) / (n + ntot);
			ntot += n;
		}
	}
	assert(ntot == recv_loglik_count(ll));
	
	return dev_mean;
}

double recv_loglik_last_dev(const struct recv_loglik *ll)
{
	assert(ll);
	assert(ll->last);
	
	return recv_sloglik_last_dev(ll->last);
}
