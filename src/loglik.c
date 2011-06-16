#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "sloglik.h"
#include "loglik.h"

static void iproc_loglik_free(iproc_loglik * loglik)
{
	if (loglik) {
		struct array *array = &loglik->sloglik_array;
		ssize_t i, n = array_count(array);

		for (i = 0; i < n; i++) {
			iproc_sloglik *sll =
			    *(iproc_sloglik **) array_item(array, i);
			iproc_sloglik_unref(sll);
		}

		array_deinit(array);
		vector_free(loglik->grad);
		model_free(loglik->model);
		free(loglik);
	}
}

static iproc_loglik *iproc_loglik_new_empty(struct model *model)
{
	iproc_loglik *loglik = calloc(1, sizeof(*loglik));
	struct design *design = model_design(model);
	ssize_t nsender = design_send_count(design);

	if (!loglik)
		return NULL;

	loglik->model = model_ref(model);
	loglik->grad = vector_alloc(design_recv_dim(design));
	loglik->grad_cached = false;
	loglik->nsend = 0;
	loglik->nrecv = 0;
	refcount_init(&loglik->refcount);
	array_init(&loglik->sloglik_array, sizeof(iproc_sloglik *));
	if (!(loglik->grad)) {
		iproc_loglik_free(loglik);
		loglik = NULL;
	}

	array_add_range(&loglik->sloglik_array, NULL, nsender);

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

static iproc_sloglik *iproc_loglik_sloglik(iproc_loglik * loglik, ssize_t isend)
{
	struct array *array = &loglik->sloglik_array;
	iproc_sloglik *sll = *(iproc_sloglik **) array_item(array, isend);

	if (!sll) {
		struct model *model = loglik->model;
		sll = iproc_sloglik_new(model, isend);
		*(iproc_sloglik **) array_item(array, isend) = sll;
	}

	return sll;
}

void iproc_loglik_insert(iproc_loglik * loglik,
			 const struct frame *f, const struct message *msg)
{
	iproc_sloglik *sll = iproc_loglik_sloglik(loglik, msg->from);
	iproc_sloglik_insertm(sll, f, msg->to, msg->nto);
	loglik->grad_cached = false;
	loglik->nsend += 1;
	loglik->nrecv += msg->nto;
}

double iproc_loglik_value(iproc_loglik * loglik)
{
	if (!loglik)
		return 0.0;

	struct array *array = &loglik->sloglik_array;
	ssize_t i, n = array_count(array);
	iproc_sloglik *sll;
	double value = 0.0;

	for (i = 0; i < n; i++) {
		sll = *(iproc_sloglik **) array_item(array, i);
		value += iproc_sloglik_value(sll);
	}

	return value;
}

static void
iproc_vector_acc_loglik_grad_nocache(struct vector *dst_vector,
				     double scale, iproc_loglik * loglik)
{
	struct array *array = &loglik->sloglik_array;
	ssize_t nsend = array_count(array);
	ssize_t i;
	iproc_sloglik *sll;

	for (i = 0; i < nsend; i++) {
		sll = *(iproc_sloglik **) array_item(array, i);
		if (sll) {
			struct vector *g = iproc_sloglik_grad(sll);
			vector_axpy(scale, g, dst_vector);
		}
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
