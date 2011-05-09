#ifndef _EVENT_TRACE_H
#define _EVENT_TRACE_H

#include "darray.h"

/* An event trace is like a mini-history.  The two main uses for a trace are
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

struct event_meta {
	double time;
	intptr_t attr;
};

struct event {
	ssize_t e;
	struct event_meta meta;
};

struct events {
	ssize_t e;
	struct darray meta;
};

struct event_trace {
	double tcur;
	struct darray pending;
	struct darray events;
};

bool event_trace_init(struct event_trace *trace);
void event_trace_deinit(struct event_trace *trace);

void event_trace_clear(struct event_trace *trace);

bool event_trace_insert(struct event_trace *trace, ssize_t e, intptr_t attr);
bool event_trace_advance_to(struct event_trace *trace, double time);
bool event_trace_reserve_insert(struct event_trace *trace, ssize_t ninsert);

double event_trace_tcur(const struct event_trace *trace);
ssize_t event_trace_size(const struct event_trace *trace);
struct events *event_trace_at(const struct event_trace *trace, ssize_t i);
struct events *event_trace_lookup(const struct event_trace *trace, ssize_t e);

ssize_t events_id(const struct events *events);
bool events_empty(const struct events *events);
ssize_t events_size(const struct events *events);
struct event_meta *events_at(const struct events *events, ssize_t i);
struct event_meta *events_back(const struct events *events);

#endif /* _EVENT_TRACE_H */
