#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>
#include <string.h>
#include <assert.h>
#include "compare.h"
#include "trace.h"
#include "memory.h"


static void
iproc_trace_clear_pending (iproc_trace *trace)
{
    assert(trace);
    iproc_array_set_size(trace->pending, 0);
}

static int64_t
iproc_trace_npending (iproc_trace *trace)
{
    assert(trace);
    return iproc_array_size(trace->pending);
}

static iproc_event *
iproc_trace_get_pending (iproc_trace *trace,
                         int64_t      i)
{
    assert(trace);
    assert(0 <= i);
    assert(i < iproc_trace_npending(trace));
    
    return &iproc_array_index(trace->pending, iproc_event, i);
}

static void
iproc_events_init (iproc_events *events,
                   int64_t       e)
{
    assert(events);
    events->e = e;
    events->meta = iproc_array_new(sizeof(iproc_event_meta));
}

static void
iproc_events_deinit (iproc_events *events)
{
    if (!events)
        return;

    iproc_array_unref(events->meta);
}

static void
iproc_events_append (iproc_events     *events,
                     iproc_event_meta *meta)
{
    assert(events);
    assert(meta);
    iproc_array_append(events->meta, meta);
}

static void
iproc_events_array_deinit (iproc_array *events_array)
{
    if (!events_array)
        return;
    
    int64_t i, n = iproc_array_size(events_array);
    for (i = 0; i < n; i++) {
        iproc_events *events = &iproc_array_index(events_array,
                                                  iproc_events,
                                                  i);
        iproc_events_deinit(events);
    }
}

static int64_t
iproc_trace_find_index (iproc_trace *trace,
                        int64_t      e)
{
    assert(trace);
    return iproc_array_bsearch(trace->events, &e, iproc_int64_compare);
}

static void
iproc_trace_free (iproc_trace *trace)
{
    if (trace) {
        iproc_events_array_deinit(trace->events);
        iproc_array_unref(trace->events);
        iproc_array_unref(trace->pending);
        iproc_free(trace);
    }
}

iproc_trace *
iproc_trace_new ()
{
    iproc_trace *trace = iproc_malloc(sizeof(*trace));

    if (!trace) return NULL;

    trace->tcur = -INFINITY;
    trace->pending = iproc_array_new(sizeof(iproc_event));
    trace->events = iproc_array_new(sizeof(iproc_events));
    iproc_refcount_init(&trace->refcount);

    if (!(trace->pending && trace->events)) {
        iproc_trace_free(trace);
        trace = NULL;
    }

    return trace;
}

iproc_trace *
iproc_trace_ref (iproc_trace *trace)
{
    if (trace) {
        iproc_refcount_get(&trace->refcount);
    }
    return trace;
}

static void
iproc_trace_release (iproc_refcount *refcount)
{
    iproc_trace *trace = container_of(refcount, iproc_trace, refcount);
    iproc_trace_free(trace);
}

void
iproc_trace_unref (iproc_trace *trace)
{
    if (!trace)
        return;

    iproc_refcount_put(&trace->refcount, iproc_trace_release);
}

void
iproc_trace_clear (iproc_trace *trace)
{
    trace->tcur = -INFINITY;
    iproc_trace_clear_pending(trace);
    iproc_array_set_size(trace->events, 0);
}

void
iproc_trace_insert (iproc_trace *trace,
                     int64_t       e)
{
    assert(trace);

    iproc_event event;
    
    event.e = e;
    event.meta.time = iproc_trace_tcur(trace);
    event.meta.attr = IPROC_EVENT_ATTR_MISSING;

    iproc_array_append(trace->pending, &event);
}

void
iproc_trace_advance_to (iproc_trace *trace,
                         double      t)
{
    assert(trace);
    assert(t >= iproc_trace_tcur(trace));

    double t0 = trace->tcur;
    
    if (t == t0)
        return;

    iproc_array *events_array = trace->events;
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
            iproc_array_insert(events_array, pos, &new_events);
        }
        
        // Add the new meta-data
        iproc_events *events = &iproc_array_index(events_array,
                                                  iproc_events,
                                                  pos);
        
        iproc_events_append(events, &event->meta);
    }

    iproc_trace_clear_pending(trace);
    trace->tcur = t;
}

double
iproc_trace_tcur (iproc_trace *trace)
{
    assert(trace);
    return trace->tcur;
}

int64_t
iproc_trace_size (iproc_trace *trace)
{
    assert(trace);
    return iproc_array_size(trace->events);
}

iproc_events *
iproc_trace_lookup (iproc_trace *trace,
                    int64_t      e)
{
    assert(trace);
    iproc_events *events = NULL;
    int64_t i = iproc_trace_find_index(trace, e);

    if (i >= 0) {
        events = &iproc_array_index(trace->events, iproc_events, i);
    }
    
    return events;
}

iproc_events *
iproc_trace_get (iproc_trace *trace,
                 int64_t      i)
{
    assert(trace);
    assert(0 <= i);
    assert(i < iproc_trace_size(trace));

    return &iproc_array_index(trace->events, iproc_events, i);
}

int64_t
iproc_events_id (iproc_events *events)
{
    assert(events);
    return events->e;
}

int64_t
iproc_events_size (iproc_events *events)
{
    assert(events);
    return iproc_array_size(events->meta);
}

iproc_event_meta *
iproc_events_get (iproc_events *events,
                  int64_t       i)
{
    assert(events);
    assert(i >= 0);
    assert(i <= iproc_events_size(events));
    
    iproc_event_meta *meta = &iproc_array_index(events->meta,
                                                iproc_event_meta,
                                                i);
    return meta;
}

iproc_event_meta *
iproc_events_last (iproc_events *events)
{
    assert(events);
    iproc_event_meta *meta = NULL;
    int64_t n = iproc_events_size(events);
    if (n > 0)
        meta = iproc_events_get(events, n - 1);
    
    return meta;
}