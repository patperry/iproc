#ifndef _IPROC_TRACE_H
#define _IPROC_TRACE_H

#include <stddef.h>
#include <stdint.h>
#include "darray.h"
#include "refcount.h"

/* A trace is like a mini-history.  The two main uses for a trace are
 * for storing the events related to each sender, and storing the events
 * related to each receiver.  For example, in the message stream
 *
 *         { (t0,1,2), (t1,1,3), (t2,3,2), (t3,1,2) },
 *
 * at time t4 > t3, the following characterizations hold:
 *
 *         sender 1's trace is { (t0,2), (t1,3), (t3,2) };
 *         sender 2's trace is { };
 *         sender 3's trace is { (t2,2) };
 *         receiver 1's trace is { };
 *         receiver 2's trace is { (t0,1), (t2,3), (t3,1) };
 *         receiver 3's trace is { (t1,1) }.
 * 
 */


#define IPROC_EVENT_ATTR_MISSING INT64_MIN
#define IPROC_EVENT_ATTR_MIN     (INT64_MIN + 1)
#define IPROC_EVENT_ATTR_MAX     INT64_MAX


typedef struct _iproc_event      iproc_event;
typedef struct _iproc_event_meta iproc_event_meta;
typedef struct _iproc_events     iproc_events;
typedef struct _iproc_trace      iproc_trace;

struct _iproc_event_meta {
    double  time;
    int64_t attr;
};

struct _iproc_event {
    int64_t          e;
    iproc_event_meta meta;
};

struct _iproc_events {
    int64_t       e;
    struct darray meta;
};

struct _iproc_trace {
    double          tcur;
    struct darray   pending;
    struct darray   events;
    struct refcount refcount;
};

iproc_trace *      iproc_trace_new        ();
iproc_trace *      iproc_trace_ref        (iproc_trace *trace);
void               iproc_trace_unref      (iproc_trace *trace);


void               iproc_trace_clear      (iproc_trace *trace);

void               iproc_trace_insert     (iproc_trace *trace,
                                           int64_t      e);
void               iproc_trace_advance_to (iproc_trace *trace,
                                           double       time);

double             iproc_trace_tcur       (iproc_trace *trace);
int64_t            iproc_trace_size       (iproc_trace *trace);
iproc_events *     iproc_trace_get        (iproc_trace *trace,
                                           int64_t      i);
iproc_events *     iproc_trace_lookup     (iproc_trace *trace,
                                           int64_t      e);

int64_t            iproc_events_id        (iproc_events *events);
int64_t            iproc_events_size      (iproc_events *events);
iproc_event_meta * iproc_events_get       (iproc_events *events,
                                           int64_t       i);
iproc_event_meta * iproc_events_last      (iproc_events *events);


#endif /* _IPROC_TRACE_H */
