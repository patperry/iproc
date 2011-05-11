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
	return list_init(&events->meta, sizeof(struct event_meta));
}

static void events_deinit(struct events *events)
{
	assert(events);
	list_deinit(&events->meta);
}

static void events_array_deinit(struct list *events_array)
{
	if (!events_array)
		return;

	ssize_t i, n = list_count(events_array);
	for (i = 0; i < n; i++) {
		struct events *events = list_item(events_array, i);
		events_deinit(events);
	}

	list_deinit(events_array);
}

static ssize_t event_trace_find_index(const struct event_trace *trace,
				      ssize_t e)
{
	assert(trace);
	return list_binary_search(&trace->events, &e, ssize_compare);
}

void event_trace_deinit(struct event_trace *trace)
{
	assert(trace);
	events_array_deinit(&trace->events);
	list_deinit(&trace->pending);
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

	if (!list_init(&trace->pending, sizeof(struct event)))
		goto fail_pending;

	if (!list_init(&trace->events, sizeof(struct events)))
		goto fail_events;

	trace->tcur = -INFINITY;
	return true;

fail_events:
	list_deinit(&trace->pending);
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
	list_clear(&trace->pending);

	n = list_count(&trace->events);
	for (i = 0; i < n; i++) {
		events_deinit(list_item(&trace->events, i));
	}

	list_clear(&trace->events);
}

bool event_trace_insert(struct event_trace *trace, ssize_t e, intptr_t attr)
{
	assert(trace);

	struct event event;

	event.e = e;
	event.meta.time = event_trace_tcur(trace);
	event.meta.attr = attr;

	return list_add(&trace->pending, &event);
}

bool event_trace_reserve_insert(struct event_trace *trace, ssize_t ninsert)
{
	assert(trace);
	assert(ninsert >= 0);
	return list_set_capacity(&trace->pending,
			      list_count(&trace->pending) + ninsert);
}

bool event_trace_advance_to(struct event_trace *trace, double t)
{
	assert(trace);
	assert(t >= event_trace_tcur(trace));

	double t0 = trace->tcur;

	if (t == t0)
		return true;

	ssize_t i, n = list_count(&trace->pending);

	// Process all pending events
	for (i = 0; i < n; i++) {
		struct event *event = list_item(&trace->pending, i);
		ssize_t pos = event_trace_find_index(trace, event->e);

		// If event doesn't already exist in events array, insert it
		if (pos < 0) {
			pos = ~pos;

			struct events *new_events;
			if ((new_events = list_insert(&trace->events, pos, NULL))) {
				if (!events_init(new_events, event->e)) {
					list_remove_at(&trace->events, pos);
					goto rollback;
				}
			}
		}
		// Add the new meta-data
		struct events *events = list_item(&trace->events, pos);

		if (!list_add(&events->meta, &event->meta))
			goto rollback;
	}
	list_clear(&trace->pending);
	trace->tcur = t;
	return true;

rollback:
	for (; i > 0; i--) {
		struct event *event = list_item(&trace->pending, i - 1);
		ssize_t pos = event_trace_find_index(trace, event->e);	// always exists
		struct events *events = list_item(&trace->events, pos);
		list_remove_at(&events->meta, list_count(&events->meta) - 1);
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
	return list_count(&trace->events);
}

struct events *event_trace_lookup(const struct event_trace *trace, ssize_t e)
{
	assert(trace);
	struct events *events = NULL;
	ssize_t i = event_trace_find_index(trace, e);

	if (i >= 0) {
		events = list_item(&trace->events, i);
	}

	return events;
}

struct events *event_trace_at(const struct event_trace *trace, ssize_t i)
{
	assert(trace);
	assert(0 <= i);
	assert(i < event_trace_size(trace));

	return list_item(&trace->events, i);
}

ssize_t events_id(const struct events *events)
{
	assert(events);
	return events->e;
}

bool events_empty(const struct events *events)
{
	assert(events);
	return !list_count(&events->meta);
}

ssize_t events_size(const struct events *events)
{
	assert(events);
	return list_count(&events->meta);
}

struct event_meta *events_at(const struct events *events, ssize_t i)
{
	assert(events);
	assert(i >= 0);
	assert(i <= events_size(events));

	return list_item(&events->meta, i);
}

struct event_meta *events_back(const struct events *events)
{
	assert(events);
	assert(!events_empty(events));
	return list_item(&events->meta, list_count(&events->meta) - 1);
}
