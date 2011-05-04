#ifndef _IPROC_HISTORY_H
#define _IPROC_HISTORY_H

#include "darray.h"
#include "event-trace.h"
#include "refcount.h"

/* A `history` object stores the history of the process up to the
 * current time.  Internally, the history maintains two arrays, one
 * containing the sender traces, and one containing the receiver traces.
 * The message (t,i,j) gets recorded as two events:
 *
 *         (1) event (t,j) in sender i's trace; and
 *         (2) event (t,i) in receiver j's trace.
 *
 * Events can be added to the history by `iproc_history_insert` and
 * `iproc_history_insertm`, the latter of which is for multicast events.
 * The current time can be changed by calling `iproc_history_advance_to`.
 *
 * The function `iproc_history_send(h,i)` gets sender i's trace;
 * the function `iproc_history_recv(h,j)` gets receiver j's trace.
 */

typedef struct _iproc_history iproc_history;
typedef struct _iproc_history_trace iproc_history_trace;

struct _iproc_history_trace {
	struct event_trace *trace;
	double tcur;
};

struct _iproc_history {
	double tcur;
	struct darray send;
	struct darray recv;
	struct refcount refcount;
};

iproc_history *iproc_history_new();
iproc_history *iproc_history_ref(iproc_history * history);
void iproc_history_unref(iproc_history * history);
void iproc_history_clear(iproc_history * history);

double iproc_history_tcur(iproc_history * history);
void iproc_history_advance_to(iproc_history * history, double t);
void iproc_history_insert(iproc_history * history, ssize_t from, ssize_t to);
void iproc_history_insertm(iproc_history * history,
			   ssize_t from, ssize_t *to, ssize_t nto);

ssize_t iproc_history_nsend(iproc_history * history);
ssize_t iproc_history_nrecv(iproc_history * history);
struct event_trace *iproc_history_send(iproc_history * history, ssize_t i);
struct event_trace *iproc_history_recv(iproc_history * history, ssize_t j);

#endif /* _IPROC_HISTORY_H */
