#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "sloglik.h"
#include "loglik.h"

static void iproc_loglik_free(iproc_loglik * loglik)
{
	if (loglik) {
		struct recv_sloglik *sll;
		ARRAY_FOREACH(sll, &loglik->slogliks) {
			recv_sloglik_deinit(sll);
		}
		array_deinit(&loglik->slogliks);
		vector_free(loglik->grad);
		model_free(loglik->model);
		free(loglik);
	}
}

static iproc_loglik *iproc_loglik_new_empty(struct model *model)
{
	iproc_loglik *loglik = calloc(1, sizeof(*loglik));
	struct design *design = model_design(model);
	ssize_t isend, nsender = design_send_count(design);

	if (!loglik)
		return NULL;

	loglik->model = model_ref(model);
	loglik->grad = vector_alloc(design_recv_dim(design));
	loglik->grad_cached = false;
	loglik->nsend = 0;
	loglik->nrecv = 0;
	refcount_init(&loglik->refcount);
	array_init(&loglik->slogliks, sizeof(struct recv_sloglik));
	
	array_set_capacity(&loglik->slogliks, nsender);
	for (isend = 0; isend < nsender; isend++) {
		struct recv_sloglik *sll = array_add(&loglik->slogliks, NULL);
		recv_sloglik_init(sll, model, isend);
	}
	
	return loglik;
}

iproc_loglik *iproc_loglik_new(struct model * model, struct messages * messages)
{
	assert(model);
	assert(!messages
	       || messages_max_from(messages) < model_sender_count(model));
	assert(!messages
	       || messages_max_to(messages) < model_receiver_count(model));

	iproc_loglik *loglik = iproc_loglik_new_empty(model);
	struct frame frame;
	frame_init(&frame, model->design);

	if (!messages)
		return loglik;

	struct messages_iter it;

	MESSAGES_FOREACH(it, messages) {
		ssize_t tie, ntie = MESSAGES_COUNT(it);

		for (tie = 0; tie < ntie; tie++) {
			const struct message *msg = MESSAGES_VAL(it, tie);
			iproc_loglik_insert(loglik, &frame, msg);
			frame_add(&frame, msg);
		}
	}

	frame_deinit(&frame);
	return loglik;
}

iproc_loglik *iproc_loglik_ref(iproc_loglik * loglik)
{
	if (loglik) {
		refcount_get(&loglik->refcount);
	}

	return loglik;
}

static void iproc_loglik_release(struct refcount *refcount)
{
	iproc_loglik *loglik = container_of(refcount, iproc_loglik, refcount);
	iproc_loglik_free(loglik);
}

void iproc_loglik_unref(iproc_loglik * loglik)
{
	if (loglik) {
		refcount_put(&loglik->refcount, iproc_loglik_release);
	}
}

static struct recv_sloglik *iproc_loglik_sloglik(iproc_loglik * loglik, ssize_t isend)
{
	return array_item(&loglik->slogliks, isend);
}

void iproc_loglik_insert(iproc_loglik * loglik,
			 const struct frame *f, const struct message *msg)
{
	struct recv_sloglik *sll = iproc_loglik_sloglik(loglik, msg->from);
	recv_sloglik_add(sll, f, msg->to, msg->nto);
	loglik->grad_cached = false;
	loglik->nsend += 1;
	loglik->nrecv += msg->nto;
}

double iproc_loglik_value(iproc_loglik * loglik)
{
	assert(loglik);
	
	double value = 0.0;
	struct recv_sloglik *sll;

	ARRAY_FOREACH(sll, &loglik->slogliks) {
		value += recv_sloglik_value(sll);
	}

	return value;
}

static void
iproc_vector_acc_loglik_grad_nocache(struct vector *dst_vector,
				     double scale, iproc_loglik * loglik)
{
	struct recv_sloglik *sll;
	
	ARRAY_FOREACH(sll, &loglik->slogliks) {
		recv_sloglik_axpy_grad(scale, sll, dst_vector);
	}
}

struct vector *iproc_loglik_grad(iproc_loglik * loglik)
{
	if (!loglik)
		return NULL;

	if (!loglik->grad_cached) {
		vector_fill(loglik->grad, 0.0);
		iproc_vector_acc_loglik_grad_nocache(loglik->grad, 1.0, loglik);
		loglik->grad_cached = true;
	}
	return loglik->grad;
}
