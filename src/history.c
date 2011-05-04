#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "history.h"

static bool trace_array_grow(struct darray *array, ssize_t n)
{
	assert(array);
	ssize_t nold = darray_size(array);
	ssize_t i;
	iproc_history_trace *ht;

	if (n > nold) {
		if (darray_resize(array, n)) {
			for (i = nold; i < n; i++) {
				ht = darray_at(array, i);
				if (!event_trace_init(&ht->trace)) {
					darray_resize(array, i);
					return false;
				}
			}
		}
	}

	return true;
}

static void trace_array_clear(struct darray *array)
{
	assert(array);
	ssize_t n = darray_size(array);
	ssize_t i;
	iproc_history_trace *ht;

	for (i = 0; i < n; i++) {
		ht = darray_at(array, i);
		ht->tcur = -INFINITY;
		event_trace_clear(&ht->trace);
	}
}

static void trace_array_deinit(struct darray *array)
{
	ssize_t n = darray_size(array);
	ssize_t i;

	for (i = 0; i < n; i++) {
		iproc_history_trace *ht = darray_at(array, i);
		event_trace_deinit(&ht->trace);
	}

	darray_deinit(array);
}

static struct event_trace *trace_array_get(double tcur,
					   struct darray *array, ssize_t i)
{
	assert(array);
	assert(i >= 0);

	if (trace_array_grow(array, i + 1)) {
		iproc_history_trace *ht = darray_at(array, i);
		struct event_trace *t = &ht->trace;

		if (ht->tcur != tcur) {
			if (!event_trace_advance_to(t, tcur))
				return NULL;
			ht->tcur = tcur;
		}
		return t;
	}
	return NULL;
}

static void iproc_history_free(iproc_history * history)
{
	if (history) {
		trace_array_deinit(&history->send);
		trace_array_deinit(&history->recv);
		free(history);
	}
}

iproc_history *iproc_history_new()
{
	iproc_history *history = calloc(1, sizeof(*history));

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

bool iproc_history_advance_to(iproc_history * history, double t)
{
	assert(history);
	assert(history->tcur <= t);

	history->tcur = t;
	return true;
}

bool iproc_history_insert(iproc_history * history, ssize_t from, ssize_t to, intptr_t attr)
{
	assert(history);
	assert(from >= 0);
	assert(to >= 0);
	return iproc_history_insertm(history, from, &to, 1, attr);
}

bool iproc_history_insertm(iproc_history * history, ssize_t from, ssize_t *to, ssize_t nto, intptr_t attr)
{
	assert(history);
	assert(to || nto == 0);
	assert(nto >= 0);

	struct event_trace *efrom = iproc_history_send(history, from);
	if (!(efrom && event_trace_reserve_insert(efrom, nto)))
		return false;
	
	ssize_t i;
	for (i = 0; i < nto; i++) {
		struct event_trace *eto = iproc_history_recv(history, to[i]);
		if (!(eto && event_trace_reserve_insert(eto, 1)))
			return false;
	}

	for (i = 0; i < nto; i++) {
		struct event_trace *eto = iproc_history_recv(history, to[i]);
		event_trace_insert(efrom, to[i], attr);
		event_trace_insert(eto, from, attr);
	}
	
	return true;
}

ssize_t iproc_history_nsend(iproc_history * history)
{
	assert(history);
	return darray_size(&history->send);
}

ssize_t iproc_history_nrecv(iproc_history * history)
{
	assert(history);
	return darray_size(&history->recv);
}

struct event_trace *iproc_history_send(iproc_history * history, ssize_t i)
{
	assert(history);
	assert(0 <= i);

	return trace_array_get(history->tcur, &history->send, i);
}

struct event_trace *iproc_history_recv(iproc_history * history, ssize_t j)
{
	assert(history);
	assert(0 <= j);

	return trace_array_get(history->tcur, &history->recv, j);
}
