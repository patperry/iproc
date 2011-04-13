#include "port.h"

#include <math.h>
#include <string.h>
#include <assert.h>
#include "compare.h"
#include "trace.h"
#include "memory.h"

static void iproc_trace_clear_pending(iproc_trace * trace)
{
	assert(trace);
	darray_resize(&trace->pending, 0);
}

static int64_t iproc_trace_npending(iproc_trace * trace)
{
	assert(trace);
	return darray_size(&trace->pending);
}

static iproc_event *iproc_trace_get_pending(iproc_trace * trace, int64_t i)
{
	assert(trace);
	assert(0 <= i);
	assert(i < iproc_trace_npending(trace));

	return &darray_index(&trace->pending, iproc_event, i);
}

static void iproc_events_init(iproc_events * events, int64_t e)
{
	assert(events);
	events->e = e;
	darray_init(&events->meta, iproc_event_meta);
}

static void iproc_events_deinit(iproc_events * events)
{
	if (!events)
		return;

	darray_deinit(&events->meta);
}

static void iproc_events_append(iproc_events * events, iproc_event_meta * meta)
{
	assert(events);
	assert(meta);
	darray_push_back(&events->meta, meta);
}

static void iproc_events_array_deinit(struct darray *events_array)
{
	if (!events_array)
		return;

	int64_t i, n = darray_size(events_array);
	for (i = 0; i < n; i++) {
		iproc_events *events = &darray_index(events_array,
						     iproc_events,
						     i);
		iproc_events_deinit(events);
	}

	darray_deinit(events_array);
}

static int64_t iproc_trace_find_index(iproc_trace * trace, int64_t e)
{
	assert(trace);
	return darray_binary_search(&trace->events, &e, int64_compare);
}

static void iproc_trace_free(iproc_trace * trace)
{
	if (trace) {
		iproc_events_array_deinit(&trace->events);
		darray_deinit(&trace->pending);
		iproc_free(trace);
	}
}

iproc_trace *iproc_trace_new()
{
	iproc_trace *trace = iproc_malloc(sizeof(*trace));

	if (!trace)
		return NULL;

	trace->tcur = -INFINITY;
	refcount_init(&trace->refcount);

	if (!(darray_init(&trace->pending, iproc_event)
	      && darray_init(&trace->events, iproc_events))) {
		iproc_trace_free(trace);
		trace = NULL;
	}

	return trace;
}

iproc_trace *iproc_trace_ref(iproc_trace * trace)
{
	if (trace) {
		refcount_get(&trace->refcount);
	}
	return trace;
}

static void iproc_trace_release(struct refcount *refcount)
{
	iproc_trace *trace = container_of(refcount, iproc_trace, refcount);
	iproc_trace_free(trace);
}

void iproc_trace_unref(iproc_trace * trace)
{
	if (!trace)
		return;

	refcount_put(&trace->refcount, iproc_trace_release);
}

void iproc_trace_clear(iproc_trace * trace)
{
	ssize_t i, n;

	trace->tcur = -INFINITY;
	iproc_trace_clear_pending(trace);

	n = darray_size(&trace->events);
	for (i = 0; i < n; i++) {
		iproc_events_deinit(&darray_index
				    (&trace->events, iproc_events, i));
	}

	darray_clear(&trace->events);
}

void iproc_trace_insert(iproc_trace * trace, int64_t e)
{
	assert(trace);

	iproc_event event;

	event.e = e;
	event.meta.time = iproc_trace_tcur(trace);
	event.meta.attr = IPROC_EVENT_ATTR_MISSING;

	darray_push_back(&trace->pending, &event);
}

void iproc_trace_advance_to(iproc_trace * trace, double t)
{
	assert(trace);
	assert(t >= iproc_trace_tcur(trace));

	double t0 = trace->tcur;

	if (t == t0)
		return;

	struct darray *events_array = &trace->events;
	int64_t i, n = iproc_trace_npending(trace);

	// Process all pending events
	for (i = 0; i < n; i++) {
		iproc_event *event = iproc_trace_get_pending(trace, i);
		int64_t pos = iproc_trace_find_index(trace, event->e);

		// If event doesn't already exist in events array, insert it
		if (pos < 0) {
			pos = ~pos;

			iproc_events new_events;
			iproc_events_init(&new_events, event->e);
			darray_insert(events_array, pos, &new_events);
		}
		// Add the new meta-data
		iproc_events *events = &darray_index(events_array,
						     iproc_events,
						     pos);

		iproc_events_append(events, &event->meta);
	}

	iproc_trace_clear_pending(trace);
	trace->tcur = t;
}

double iproc_trace_tcur(iproc_trace * trace)
{
	assert(trace);
	return trace->tcur;
}

int64_t iproc_trace_size(iproc_trace * trace)
{
	assert(trace);
	return darray_size(&trace->events);
}

iproc_events *iproc_trace_lookup(iproc_trace * trace, int64_t e)
{
	assert(trace);
	iproc_events *events = NULL;
	int64_t i = iproc_trace_find_index(trace, e);

	if (i >= 0) {
		events = &darray_index(&trace->events, iproc_events, i);
	}

	return events;
}

iproc_events *iproc_trace_get(iproc_trace * trace, int64_t i)
{
	assert(trace);
	assert(0 <= i);
	assert(i < iproc_trace_size(trace));

	return &darray_index(&trace->events, iproc_events, i);
}

int64_t iproc_events_id(iproc_events * events)
{
	assert(events);
	return events->e;
}

int64_t iproc_events_size(iproc_events * events)
{
	assert(events);
	return darray_size(&events->meta);
}

iproc_event_meta *iproc_events_get(iproc_events * events, int64_t i)
{
	assert(events);
	assert(i >= 0);
	assert(i <= iproc_events_size(events));

	iproc_event_meta *meta = &darray_index(&events->meta,
					       iproc_event_meta,
					       i);
	return meta;
}

iproc_event_meta *iproc_events_last(iproc_events * events)
{
	assert(events);
	iproc_event_meta *meta = NULL;
	int64_t n = iproc_events_size(events);
	if (n > 0)
		meta = iproc_events_get(events, n - 1);

	return meta;
}
