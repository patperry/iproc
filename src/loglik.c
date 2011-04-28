#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "sloglik.h"
#include "loglik.h"

static void iproc_loglik_free(iproc_loglik * loglik)
{
	if (loglik) {
		struct darray *array = &loglik->sloglik_array;
		ssize_t i, n = darray_size(array);

		for (i = 0; i < n; i++) {
			iproc_sloglik *sll =
			    *(iproc_sloglik **) darray_at(array, i);
			iproc_sloglik_unref(sll);
		}

		darray_deinit(array);
		vector_free(loglik->grad);
		iproc_model_unref(loglik->model);
		free(loglik);
	}
}

static iproc_loglik *iproc_loglik_new_empty(iproc_model * model)
{
	iproc_loglik *loglik = calloc(1, sizeof(*loglik));
	iproc_design *design = iproc_model_design(model);
	ssize_t nsender = iproc_design_nsender(design);

	if (!loglik)
		return NULL;

	loglik->model = iproc_model_ref(model);
	loglik->grad = vector_alloc(iproc_design_dim(design));
	loglik->grad_cached = false;
	loglik->nsend = 0;
	loglik->nrecv = 0;
	refcount_init(&loglik->refcount);

	if (!(darray_init(&loglik->sloglik_array, sizeof(iproc_sloglik *))
	      && loglik->grad)) {
		iproc_loglik_free(loglik);
		loglik = NULL;
	}

	darray_resize(&loglik->sloglik_array, nsender);

	return loglik;
}

iproc_loglik *iproc_loglik_new(iproc_model * model, struct messages * messages)
{
	assert(model);
	assert(!messages
	       || messages_max_from(messages) <
	       iproc_model_nsender(model));
	assert(!messages
	       || messages_max_to(messages) <
	       iproc_model_nreceiver(model));

	iproc_loglik *loglik = iproc_loglik_new_empty(model);

	if (!messages)
		return loglik;

	struct messages_iter *it = messages_iter_alloc(messages);

	while (messages_iter_next(it)) {
		iproc_history *history = messages_iter_history(it);
		ssize_t tie, ntie = messages_iter_ntie(it);

		for (tie = 0; tie < ntie; tie++) {
			messages_iter_select(it, tie);

			ssize_t from = messages_iter_from(it);
			ssize_t *to = messages_iter_to(it);
			ssize_t nto = messages_iter_nto(it);

			iproc_loglik_insertm(loglik, history, from, to, nto);
		}
	}

	messages_iter_free(it);
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
	struct darray *array = &loglik->sloglik_array;
	iproc_sloglik *sll = *(iproc_sloglik **) darray_at(array, isend);

	if (!sll) {
		iproc_model *model = loglik->model;
		sll = iproc_sloglik_new(model, isend);
		*(iproc_sloglik **) darray_at(array, isend) = sll;
	}

	return sll;
}

void
iproc_loglik_insert(iproc_loglik * loglik,
		    iproc_history * history, ssize_t from, ssize_t to)
{
	iproc_sloglik *sll = iproc_loglik_sloglik(loglik, from);
	iproc_sloglik_insert(sll, history, to);
	loglik->nsend++;
	loglik->nrecv++;
}

void
iproc_loglik_insertm(iproc_loglik * loglik,
		     iproc_history * history,
		     ssize_t from, ssize_t *to, ssize_t nto)
{
	iproc_sloglik *sll = iproc_loglik_sloglik(loglik, from);
	iproc_sloglik_insertm(sll, history, to, nto);
	loglik->grad_cached = false;
	loglik->nsend += 1;
	loglik->nrecv += nto;
}

double iproc_loglik_value(iproc_loglik * loglik)
{
	if (!loglik)
		return 0.0;

	struct darray *array = &loglik->sloglik_array;
	ssize_t i, n = darray_size(array);
	iproc_sloglik *sll;
	double value = 0.0;

	for (i = 0; i < n; i++) {
		sll = *(iproc_sloglik **) darray_at(array, i);
		value += iproc_sloglik_value(sll);
	}

	return value;
}

static void
iproc_vector_acc_loglik_grad_nocache(struct vector *dst_vector,
				     double scale, iproc_loglik * loglik)
{
	struct darray *array = &loglik->sloglik_array;
	ssize_t nsend = darray_size(array);
	ssize_t i;
	iproc_sloglik *sll;

	for (i = 0; i < nsend; i++) {
		sll = *(iproc_sloglik **) darray_at(array, i);
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
