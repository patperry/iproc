#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "history.h"

static bool trace_array_grow(struct array *array, ssize_t n)
{
	assert(array);
	ssize_t nold = array_count(array);
	ssize_t i;
	struct history_trace *ht;

	if (n > nold) {
		if (array_set_capacity(array, n)) {
			for (i = nold; i < n; i++) {
				ht = array_insert(array, i, NULL);
				if (!event_trace_init(&ht->trace)) {
					array_remove_at(array, i);
					return false;
				}
			}
		}
	}

	return true;
}

static void trace_array_clear(struct array *array)
{
	assert(array);
	ssize_t n = array_count(array);
	ssize_t i;
	struct history_trace *ht;

	for (i = 0; i < n; i++) {
		ht = array_item(array, i);
		ht->tcur = -INFINITY;
		event_trace_clear(&ht->trace);
	}
}

static void trace_array_deinit(struct array *array)
{
	ssize_t n = array_count(array);
	ssize_t i;

	for (i = 0; i < n; i++) {
		struct history_trace *ht = array_item(array, i);
		event_trace_deinit(&ht->trace);
	}

	array_deinit(array);
}

static struct event_trace *trace_array_get(double tcur,
					   struct array *array, ssize_t i)
{
	assert(array);
	assert(i >= 0);

	if (trace_array_grow(array, i + 1)) {
		struct history_trace *ht = array_item(array, i);
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
	trace_array_deinit(&history->recv);
	trace_array_deinit(&history->send);
}

bool history_init(struct history *history)
{
	assert(history);

	if (!array_init(&history->send, sizeof(struct history_trace)))
		goto fail_send;

	if (!array_init(&history->recv, sizeof(struct history_trace)))
		goto fail_recv;

	history->tcur = -INFINITY;
	return true;

	array_deinit(&history->recv);
fail_recv:
	array_deinit(&history->send);
fail_send:
	return false;
}

void history_clear(struct history *history)
{
	assert(history);
	history->tcur = -INFINITY;
	trace_array_clear(&history->send);
	trace_array_clear(&history->recv);
}

double history_tcur(const struct history *history)
{
	assert(history);
	return history->tcur;
}

bool history_advance_to(struct history *history, double t)
{
	assert(history);
	assert(history->tcur <= t);

	history->tcur = t;
	return true;
}

bool history_insert(struct history *history, ssize_t from, ssize_t *to,
		    ssize_t nto, intptr_t attr)
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

ssize_t history_nsend(struct history *history)
{
	assert(history);
	return array_count(&history->send);
}

ssize_t history_nrecv(struct history *history)
{
	assert(history);
	return array_count(&history->recv);
}

struct event_trace *history_send(struct history *history, ssize_t i)
{
	assert(history);
	assert(0 <= i);

	return trace_array_get(history->tcur, &history->send, i);
}

struct event_trace *history_recv(struct history *history, ssize_t j)
{
	assert(history);
	assert(0 <= j);

	return trace_array_get(history->tcur, &history->recv, j);
}
