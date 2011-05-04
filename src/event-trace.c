#include "port.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "compare.h"
#include "event-trace.h"

DEFINE_COMPARE_FN(ssize_compare, ssize_t)

static bool events_init(struct events * events, ssize_t e)
{
	assert(events);
	events->e = e;
	return darray_init(&events->meta, sizeof(struct event_meta));
}

static void events_deinit(struct events * events)
{
	assert(events);
	darray_deinit(&events->meta);
}

static void events_array_deinit(struct darray *events_array)
{
	if (!events_array)
		return;

	ssize_t i, n = darray_size(events_array);
	for (i = 0; i < n; i++) {
		struct events *events = darray_at(events_array, i);
		events_deinit(events);
	}

	darray_deinit(events_array);
}

static ssize_t event_trace_find_index(const struct event_trace * trace, ssize_t e)
{
	assert(trace);
	return darray_binary_search(&trace->events, &e, ssize_compare);
}

void event_trace_deinit(struct event_trace *trace)
{
	assert(trace);
	events_array_deinit(&trace->events);
	darray_deinit(&trace->pending);
}

void event_trace_free(struct event_trace * trace)
{
	if (trace) {
		event_trace_deinit(trace);
		free(trace);
	}
}

bool event_trace_init(struct event_trace *trace)
{
	assert(trace);
	
	if (!darray_init(&trace->pending, sizeof(struct event)))
		goto fail_pending;
	
	if (!darray_init(&trace->events, sizeof(struct events)))
		goto fail_events;
	
	trace->tcur = -INFINITY;
	return true;

fail_events:
	darray_deinit(&trace->pending);
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


void event_trace_clear(struct event_trace * trace)
{
	ssize_t i, n;

	trace->tcur = -INFINITY;
	darray_clear(&trace->pending);

	n = darray_size(&trace->events);
	for (i = 0; i < n; i++) {
		events_deinit(darray_at(&trace->events, i));
	}

	darray_clear(&trace->events);
}

bool event_trace_insert(struct event_trace * trace, ssize_t e, intptr_t attr)
{
	assert(trace);

	struct event event;

	event.e = e;
	event.meta.time = event_trace_tcur(trace);
	event.meta.attr = attr;

	return darray_push_back(&trace->pending, &event);
}

bool event_trace_advance_to(struct event_trace * trace, double t)
{
	assert(trace);
	assert(t >= event_trace_tcur(trace));

	double t0 = trace->tcur;

	if (t == t0)
		return true;

	struct darray *events_array = &trace->events;
	ssize_t i, n = darray_size(&trace->pending);

	// Process all pending events
	for (i = 0; i < n; i++) {
		struct event *event = darray_at(&trace->pending, i);
		ssize_t pos = event_trace_find_index(trace, event->e);

		// If event doesn't already exist in events array, insert it
		if (pos < 0) {
			pos = ~pos;

			struct events new_events;
			if (events_init(&new_events, event->e)) {
				if (!darray_insert(events_array, pos, &new_events)) {
					events_deinit(&new_events);
					goto rollback;
				}
			}
		}
		// Add the new meta-data
		struct events *events = darray_at(events_array, pos);

		if (!darray_push_back(&events->meta, &event->meta))
			goto rollback;
	}
	darray_clear(&trace->pending);
	trace->tcur = t;
	return true;
	
rollback:
	for (; i > 0; i--) {
		struct event *event = darray_at(&trace->pending, i - 1);
		ssize_t pos = event_trace_find_index(trace, event->e); // always exists
		struct events *events = darray_at(events_array, pos);
		darray_pop_back(&events->meta);
	}
	return false;
}

double event_trace_tcur(const struct event_trace * trace)
{
	assert(trace);
	return trace->tcur;
}

ssize_t event_trace_size(const struct event_trace * trace)
{
	assert(trace);
	return darray_size(&trace->events);
}

struct events *event_trace_lookup(const struct event_trace * trace, ssize_t e)
{
	assert(trace);
	struct events *events = NULL;
	ssize_t i = event_trace_find_index(trace, e);

	if (i >= 0) {
		events = darray_at(&trace->events, i);
	}

	return events;
}

struct events *event_trace_at(const struct event_trace * trace, ssize_t i)
{
	assert(trace);
	assert(0 <= i);
	assert(i < event_trace_size(trace));

	return darray_at(&trace->events, i);
}

ssize_t events_id(const struct events * events)
{
	assert(events);
	return events->e;
}

bool events_empty(const struct events *events)
{
	assert(events);
	return darray_empty(&events->meta);
}

ssize_t events_size(const struct events * events)
{
	assert(events);
	return darray_size(&events->meta);
}

struct event_meta *events_at(const struct events * events, ssize_t i)
{
	assert(events);
	assert(i >= 0);
	assert(i <= events_size(events));

	return darray_at(&events->meta, i);
}

struct event_meta *events_back(const struct events * events)
{
	assert(events);
	assert(!events_empty(events));
	return darray_back(&events->meta);
}
