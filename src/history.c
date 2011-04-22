#include "port.h"

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "history.h"

static void trace_array_grow(struct darray *array, int64_t n)
{
	assert(array);
	int64_t nold = darray_size(array);

	if (n > nold) {
		darray_resize(array, n);
	}
}

static void trace_array_clear(struct darray *array)
{
	assert(array);
	int64_t n = darray_size(array);
	int64_t i;
	iproc_history_trace *ht;

	for (i = 0; i < n; i++) {
		ht = darray_at(array, i);
		ht->tcur = -INFINITY;

		if (ht->trace)
			iproc_trace_clear(ht->trace);
	}
}

static void trace_array_deinit(struct darray *array)
{
	ssize_t n = darray_size(array);
	ssize_t i;

	for (i = 0; i < n; i++) {
		iproc_history_trace *ht =
		    darray_at(array, i);
		if (ht->trace)
			iproc_trace_unref(ht->trace);
	}

	darray_deinit(array);
}

static iproc_trace *trace_array_get(double tcur,
				    struct darray *array, int64_t i)
{
	assert(array);
	assert(i >= 0);

	trace_array_grow(array, i + 1);

	iproc_history_trace *ht = darray_at(array,
						i);
	iproc_trace *t;

	if (!(ht->trace)) {
		ht->trace = iproc_trace_new(tcur);
	}

	t = ht->trace;

	if (ht->tcur != tcur) {
		iproc_trace_advance_to(t, tcur);
		ht->tcur = tcur;
	}

	return t;
}

static void iproc_history_free(iproc_history * history)
{
	if (history) {
		trace_array_deinit(&history->send);
		trace_array_deinit(&history->recv);
		iproc_free(history);
	}
}

iproc_history *iproc_history_new()
{
	iproc_history *history = iproc_calloc(1, sizeof(*history));

	if (history && darray_init(&history->send, sizeof(iproc_history_trace))
	    && darray_init(&history->recv, sizeof(iproc_history_trace))
	    && refcount_init(&history->refcount)) {
		history->tcur = -INFINITY;
		return history;
	}

	iproc_history_free(history);
	return NULL;
}

iproc_history *iproc_history_ref(iproc_history * history)
{
	if (history) {
		refcount_get(&history->refcount);
	}
	return history;
}

static void iproc_history_release(struct refcount *refcount)
{
	iproc_history *history =
	    container_of(refcount, iproc_history, refcount);
	iproc_history_free(history);
}

void iproc_history_unref(iproc_history * history)
{
	if (!history)
		return;

	refcount_put(&history->refcount, iproc_history_release);
}

void iproc_history_clear(iproc_history * history)
{
	assert(history);
	history->tcur = -INFINITY;
	trace_array_clear(&history->send);
	trace_array_clear(&history->recv);
}

double iproc_history_tcur(iproc_history * history)
{
	assert(history);
	return history->tcur;
}

void iproc_history_advance_to(iproc_history * history, double t)
{
	assert(history);
	assert(history->tcur <= t);

	history->tcur = t;
}

void iproc_history_insert(iproc_history * history, int64_t from, int64_t to)
{
	assert(history);
	assert(from >= 0);
	assert(to >= 0);

	iproc_trace *efrom = iproc_history_send(history, from);
	iproc_trace *eto = iproc_history_recv(history, to);

	iproc_trace_insert(efrom, to);
	iproc_trace_insert(eto, from);
}

void
iproc_history_insertm(iproc_history * history,
		      int64_t from, int64_t * to, int64_t nto)
{
	assert(history);
	assert(to || nto == 0);
	assert(nto >= 0);

	int i;

	for (i = 0; i < nto; i++) {
		iproc_history_insert(history, from, to[i]);
	}
}

int64_t iproc_history_nsend(iproc_history * history)
{
	assert(history);
	return darray_size(&history->send);
}

int64_t iproc_history_nrecv(iproc_history * history)
{
	assert(history);
	return darray_size(&history->recv);
}

iproc_trace *iproc_history_send(iproc_history * history, int64_t i)
{
	assert(history);
	assert(0 <= i);

	return trace_array_get(history->tcur, &history->send, i);
}

iproc_trace *iproc_history_recv(iproc_history * history, int64_t j)
{
	assert(history);
	assert(0 <= j);

	return trace_array_get(history->tcur, &history->recv, j);
}
