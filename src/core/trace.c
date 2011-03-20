#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>
#include <string.h>
#include <assert.h>
#include "trace.h"
#include "memory.h"

static int
compare_int64 (void *px, void *py)
{
    int64_t x = *(int64_t *)px;
    int64_t y = *(int64_t *)py;

    if (x < y) {
        return -1;
    } else if (x > y) {
        return 1;
    } else {
        return 0;
    }
}

static void
iproc_trace_clear_cur (iproc_trace *trace)
{
    assert(trace);
    iproc_array_set_size(trace->cur, 0);
}

static int64_t
iproc_trace_ncur (iproc_trace *trace)
{
    assert(trace);
    return iproc_array_size(trace->cur);
}

static int64_t
iproc_trace_cur (iproc_trace *trace,
                  int64_t       i)
{
    assert(trace);
    assert(0 <= i);
    assert(i < iproc_trace_ncur(trace));
    
    return iproc_array_index(trace->cur, int64_t, i);
}

static int64_t
iproc_trace_find_cur (iproc_trace *trace,
                       int64_t       e)
{
    assert(trace);
    return iproc_array_lfind(trace->cur, &e, compare_int64);
}

static int64_t
iproc_trace_find_index (iproc_trace *trace,
                         int64_t       e)
{
    assert(trace);
    
    /* use that iproc_past_event inherits from iproc_event */
    return iproc_array_bsearch(trace->past, &e, compare_int64);
}

static void
iproc_trace_clear_past (iproc_trace *trace)
{
    assert(trace);
    iproc_array_set_size(trace->past, 0);
}

static void
iproc_trace_free (iproc_trace *trace)
{
    if (trace) {
        iproc_array_unref(trace->cur);
        iproc_array_unref(trace->past);
        iproc_free(trace);
    }
}

iproc_trace *
iproc_trace_new ()
{
    iproc_trace *trace = iproc_malloc(sizeof(*trace));

    if (!trace) return NULL;

    trace->tcur = -INFINITY;
    trace->cur  = iproc_array_new(sizeof(int64_t));
    trace->past = iproc_array_new(sizeof(iproc_event));
    iproc_refcount_init(&trace->refcount);

    if (!(trace->cur && trace->past)) {
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
    iproc_trace_clear_cur(trace);
    iproc_trace_clear_past(trace);
}

void
iproc_trace_insert (iproc_trace *trace,
                     int64_t       e)
{
    assert(trace);

    if (iproc_trace_find_cur(trace, e) < 0) {
        iproc_array_append(trace->cur, &e);
    }
}

void
iproc_trace_advance_to (iproc_trace *trace,
                         double        t)
{
    assert(trace);
    assert(t >= trace->tcur);

    if (t == trace->tcur)
        return;

    double t0 = trace->tcur;
    int64_t e;
    int64_t ic, ip;
    int64_t nc = iproc_trace_ncur(trace);
    iproc_array *past = trace->past;

    /* Move current trace to past event set */
    for (ic = 0; ic < nc; ic++) {
        e = iproc_trace_cur(trace, ic);
        ip = iproc_trace_find_index(trace, e);

        /* If event doesn't already exist in past set, insert it */
        if (ip < 0) {
            iproc_event past_event = { e, t0 };
            iproc_array_insert(past, ~ip, &past_event);

        /* Set the time of the event to the old time */
        } else {
            iproc_array_index(past, iproc_event, ip).t = t0;
        }
    }

    iproc_trace_clear_cur(trace);
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
    return iproc_array_size(trace->past);
}

iproc_event *
iproc_trace_lookup (iproc_trace *trace,
                     int64_t       e)
{
    assert(trace);
    iproc_event *event = NULL;
    int64_t i = iproc_trace_find_index(trace, e);

    if (i >= 0) {
        event = &iproc_array_index(trace->past, iproc_event, i);
    }
    
    return event;
}

iproc_event *
iproc_trace_get (iproc_trace *trace,
                   int64_t       i)
{
    assert(trace);
    assert(0 <= i);
    assert(i < iproc_trace_size(trace));

    return &iproc_array_index(trace->past, iproc_event, i);
}
