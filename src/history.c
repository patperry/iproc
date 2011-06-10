#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "history.h"

static void trace_array_grow(struct array *array, ssize_t n)
{
	assert(array);
	ssize_t nold = array_count(array);
	ssize_t i;
	struct history_trace *ht;

	for (i = nold; i < n; i++) {
		ht = array_insert(array, i, NULL);
		event_trace_init(&ht->trace);
	}
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

	trace_array_grow(array, i + 1);
	struct history_trace *ht = array_item(array, i);

	if (ht->tcur != tcur) {
		event_trace_advance_to(&ht->trace, tcur);
		ht->tcur = tcur;
	}
	return &ht->trace;
}

void history_deinit(struct history *history)
{
	assert(history);
	trace_array_deinit(&history->recv);
	trace_array_deinit(&history->send);
}

void history_init(struct history *history)
{
	assert(history);

	array_init(&history->send, sizeof(struct history_trace));
	array_init(&history->recv, sizeof(struct history_trace));
	history->tcur = -INFINITY;
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

void history_advance_to(struct history *history, double t)
{
	assert(history);
	assert(history->tcur <= t);

	history->tcur = t;
}

void history_insert(struct history *history, ssize_t from, ssize_t *to,
		    ssize_t nto, intptr_t attr)
{
	assert(history);
	assert(to || nto == 0);
	assert(nto >= 0);

	struct event_trace *efrom = history_send(history, from);
	ssize_t i;

	for (i = 0; i < nto; i++) {
		struct event_trace *eto = history_recv(history, to[i]);
		event_trace_insert(efrom, to[i], attr);
		event_trace_insert(eto, from, attr);
	}
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
