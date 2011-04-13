#include "port.h"

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "messages.h"

static void iproc_messages_free(iproc_messages * msgs)
{
	if (msgs) {
		darray_deinit(&msgs->array);
		darray_deinit(&msgs->recipients);
		iproc_free(msgs);
	}
}

iproc_messages *iproc_messages_new()
{
	iproc_messages *msgs = iproc_calloc(1, sizeof(*msgs));
	if (!msgs)
		return NULL;

	msgs->tcur = -INFINITY;
	msgs->max_to = -1;
	msgs->max_from = -1;
	msgs->max_nto = 0;
	refcount_init(&msgs->refcount);

	if (!(darray_init(&msgs->array, iproc_message)
	      && darray_init(&msgs->recipients, int64_t))) {
		iproc_messages_free(msgs);
		msgs = NULL;
	}

	return msgs;
}

iproc_messages *iproc_messages_ref(iproc_messages * msgs)
{
	if (msgs) {
		refcount_get(&msgs->refcount);
	}
	return msgs;
}

static void iproc_messages_release(struct refcount *refcount)
{
	iproc_messages *msgs = container_of(refcount, iproc_messages, refcount);
	iproc_messages_free(msgs);
}

void iproc_messages_unref(iproc_messages * msgs)
{
	if (msgs) {
		refcount_put(&msgs->refcount, iproc_messages_release);
	}
}

int64_t iproc_messages_size(iproc_messages * msgs)
{
	assert(msgs);
	return darray_size(&msgs->array);
}

void iproc_messages_advance_to(iproc_messages * msgs, double t)
{
	assert(msgs);
	assert(t >= msgs->tcur);
	msgs->tcur = t;
}

void iproc_messages_insert(iproc_messages * msgs, int64_t from, int64_t to)
{
	assert(msgs);
	assert(from >= 0);
	assert(to >= 0);
	iproc_messages_insertm(msgs, from, &to, 1);
}

void
iproc_messages_insertm(iproc_messages * msgs,
		       int64_t from, int64_t * to, int64_t nto)
{
	assert(msgs);
	assert(from >= 0);
	assert(nto >= 0);
	assert(to || nto == 0);

	double time = msgs->tcur;
	struct darray *array = &msgs->array;
	struct darray *recipients = &msgs->recipients;

	int64_t ito = darray_size(recipients);
	iproc_message m = { time, from, ito, nto };
	int64_t i;

	for (i = 0; i < nto; i++) {
		assert(to[i] >= 0);

		darray_push_back(recipients, to + i);
		if (to[i] > msgs->max_to)
			msgs->max_to = to[i];
	}

	if (from > msgs->max_from)
		msgs->max_from = from;
	if (nto > msgs->max_nto)
		msgs->max_nto = nto;

	darray_push_back(array, &m);
}

int64_t iproc_messages_max_from(iproc_messages * msgs)
{
	assert(msgs);
	return msgs->max_from;
}

int64_t iproc_messages_max_to(iproc_messages * msgs)
{
	assert(msgs);
	return msgs->max_to;
}

int64_t iproc_messages_max_nto(iproc_messages * msgs)
{
	assert(msgs);
	return msgs->max_nto;
}
