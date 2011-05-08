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
	struct history_trace *ht;

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
	struct history_trace *ht;

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
		struct history_trace *ht = darray_at(array, i);
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
		struct history_trace *ht = darray_at(array, i);
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

void history_deinit(struct history *history)
{
	assert(history);
	refcount_deinit(&history->refcount);
	trace_array_deinit(&history->recv);	
	trace_array_deinit(&history->send);
}

bool history_init(struct history *history)
{
	assert(history);
	
	if (!darray_init(&history->send, sizeof(struct history_trace)))
		goto fail_send;
	
	if (!darray_init(&history->recv, sizeof(struct history_trace)))
		goto fail_recv;
	
	if (!refcount_init(&history->refcount))
		goto fail_refcount;
	
	history->tcur = -INFINITY;
	return true;
	
fail_refcount:
	darray_deinit(&history->recv);
fail_recv:
	darray_deinit(&history->send);
fail_send:
	return false;
}

struct history *history_alloc(void)
{
	struct history *history = malloc(sizeof(*history));

	if (history) {
		if (history_init(history))
			return history;
		free(history);
	}
	
	return NULL;
}

struct history *history_ref(struct history * history)
{
	if (history) {
		refcount_get(&history->refcount);
	}
	return history;
}

void history_free(struct history * history)
{
	if (!history)
		return;
	
	if (refcount_put(&history->refcount, NULL)) {
		history_deinit(history);
		free(history);
	}
}

void history_clear(struct history * history)
{
	assert(history);
	history->tcur = -INFINITY;
	trace_array_clear(&history->send);
	trace_array_clear(&history->recv);
}

double history_tcur(const struct history * history)
{
	assert(history);
	return history->tcur;
}

bool history_advance_to(struct history * history, double t)
{
	assert(history);
	assert(history->tcur <= t);

	history->tcur = t;
	return true;
}

bool history_insert(struct history * history, ssize_t from, ssize_t *to, ssize_t nto, intptr_t attr)
{
	assert(history);
	assert(to || nto == 0);
	assert(nto >= 0);

	struct event_trace *efrom = history_send(history, from);
	if (!(efrom && event_trace_reserve_insert(efrom, nto)))
		return false;
	
	ssize_t i;
	for (i = 0; i < nto; i++) {
		struct event_trace *eto = history_recv(history, to[i]);
		if (!(eto && event_trace_reserve_insert(eto, 1)))
			return false;
	}

	for (i = 0; i < nto; i++) {
		struct event_trace *eto = history_recv(history, to[i]);
		event_trace_insert(efrom, to[i], attr);
		event_trace_insert(eto, from, attr);
	}
	
	return true;
}

ssize_t history_nsend(struct history * history)
{
	assert(history);
	return darray_size(&history->send);
}

ssize_t history_nrecv(struct history * history)
{
	assert(history);
	return darray_size(&history->recv);
}

struct event_trace *history_send(struct history * history, ssize_t i)
{
	assert(history);
	assert(0 <= i);

	return trace_array_get(history->tcur, &history->send, i);
}

struct event_trace *history_recv(struct history * history, ssize_t j)
{
	assert(history);
	assert(0 <= j);

	return trace_array_get(history->tcur, &history->recv, j);
}
