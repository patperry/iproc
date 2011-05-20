#include "port.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "compare.h"
#include "event-trace.h"

DEFINE_COMPARE_FN(ssize_compare, ssize_t)

static bool events_init(struct events *events, ssize_t e)
{
	assert(events);
	events->e = e;
	return array_init(&events->meta, sizeof(struct event_meta));
}

static void events_deinit(struct events *events)
{
	assert(events);
	array_deinit(&events->meta);
}

static void events_array_deinit(struct array *events_array)
{
	if (!events_array)
		return;

	ssize_t i, n = array_count(events_array);
	for (i = 0; i < n; i++) {
		struct events *events = array_item(events_array, i);
		events_deinit(events);
	}

	array_deinit(events_array);
}

static ssize_t event_trace_find_index(const struct event_trace *trace,
				      ssize_t e)
{
	assert(trace);
	return array_binary_search(&trace->events, &e, ssize_compare, NULL);
}

void event_trace_deinit(struct event_trace *trace)
{
	assert(trace);
	events_array_deinit(&trace->events);
	array_deinit(&trace->pending);
}

void event_trace_free(struct event_trace *trace)
{
	if (trace) {
		event_trace_deinit(trace);
		free(trace);
	}
}

bool event_trace_init(struct event_trace *trace)
{
	assert(trace);

	if (!array_init(&trace->pending, sizeof(struct event)))
		goto fail_pending;

	if (!array_init(&trace->events, sizeof(struct events)))
		goto fail_events;

	trace->tcur = -INFINITY;
	return true;

fail_events:
	array_deinit(&trace->pending);
fail_pending:
	return false;
}

struct event_trace *event_trace_alloc()
{
	struct event_trace *trace = malloc(sizeof(*trace));

	if (trace) {
		if (event_trace_init(trace))
			return trace;
		free(trace);
	}
	return NULL;
}

void event_trace_clear(struct event_trace *trace)
{
	ssize_t i, n;

	trace->tcur = -INFINITY;
	array_clear(&trace->pending);

	n = array_count(&trace->events);
	for (i = 0; i < n; i++) {
		events_deinit(array_item(&trace->events, i));
	}

	array_clear(&trace->events);
}

bool event_trace_insert(struct event_trace *trace, ssize_t e, intptr_t attr)
{
	assert(trace);

	struct event event;

	event.e = e;
	event.meta.time = event_trace_tcur(trace);
	event.meta.attr = attr;

	return array_add(&trace->pending, &event);
}

bool event_trace_reserve_insert(struct event_trace *trace, ssize_t ninsert)
{
	assert(trace);
	assert(ninsert >= 0);
	return array_set_capacity(&trace->pending,
				  array_count(&trace->pending) + ninsert);
}

bool event_trace_advance_to(struct event_trace *trace, double t)
{
	assert(trace);
	assert(t >= event_trace_tcur(trace));

	double t0 = trace->tcur;

	if (t == t0)
		return true;

	ssize_t i, n = array_count(&trace->pending);

	// Process all pending events
	for (i = 0; i < n; i++) {
		struct event *event = array_item(&trace->pending, i);
		ssize_t pos = event_trace_find_index(trace, event->e);

		// If event doesn't already exist in events array, insert it
		if (pos < 0) {
			pos = ~pos;

			struct events *new_events;
			if ((new_events =
			     array_insert(&trace->events, pos, NULL))) {
				if (!events_init(new_events, event->e)) {
					array_remove_at(&trace->events, pos);
					goto rollback;
				}
			}
		}
		// Add the new meta-data
		struct events *events = array_item(&trace->events, pos);

		if (!array_add(&events->meta, &event->meta))
			goto rollback;
	}
	array_clear(&trace->pending);
	trace->tcur = t;
	return true;

rollback:
	for (; i > 0; i--) {
		struct event *event = array_item(&trace->pending, i - 1);
		ssize_t pos = event_trace_find_index(trace, event->e);	// always exists
		struct events *events = array_item(&trace->events, pos);
		array_remove_at(&events->meta, array_count(&events->meta) - 1);
	}
	return false;
}

double event_trace_tcur(const struct event_trace *trace)
{
	assert(trace);
	return trace->tcur;
}

ssize_t event_trace_size(const struct event_trace *trace)
{
	assert(trace);
	return array_count(&trace->events);
}

struct events *event_trace_lookup(const struct event_trace *trace, ssize_t e)
{
	assert(trace);
	struct events *events = NULL;
	ssize_t i = event_trace_find_index(trace, e);

	if (i >= 0) {
		events = array_item(&trace->events, i);
	}

	return events;
}

struct events *event_trace_at(const struct event_trace *trace, ssize_t i)
{
	assert(trace);
	assert(0 <= i);
	assert(i < event_trace_size(trace));

	return array_item(&trace->events, i);
}

ssize_t events_id(const struct events *events)
{
	assert(events);
	return events->e;
}

bool events_empty(const struct events *events)
{
	assert(events);
	return !array_count(&events->meta);
}

ssize_t events_size(const struct events *events)
{
	assert(events);
	return array_count(&events->meta);
}

struct event_meta *events_at(const struct events *events, ssize_t i)
{
	assert(events);
	assert(i >= 0);
	assert(i <= events_size(events));

	return array_item(&events->meta, i);
}

struct event_meta *events_back(const struct events *events)
{
	assert(events);
	assert(!events_empty(events));
	return array_item(&events->meta, array_count(&events->meta) - 1);
}
