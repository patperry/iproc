#ifndef _HISTORY_H
#define _HISTORY_H

#include "array.h"
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
 * Events can be added to the history by `history_insert`.
 * The current time can be changed by calling `history_advance_to`.
 *
 * The function `history_send(h,i)` gets sender i's trace;
 * the function `history_recv(h,j)` gets receiver j's trace.
 */

struct history_trace {
	struct event_trace trace;
	double tcur;
};

struct history {
	double tcur;
	struct array send;
	struct array recv;
};

void history_init(struct history *history);
void history_deinit(struct history *history);

void history_clear(struct history *history);

double history_tcur(const struct history *history);
void history_advance_to(struct history *history, double t);
void history_add(struct history *history, ssize_t from, ssize_t *to,
		    ssize_t nto, intptr_t attr);

ssize_t history_nsend(struct history *history);
ssize_t history_nrecv(struct history *history);
struct event_trace *history_send(struct history *history, ssize_t i);
struct event_trace *history_recv(struct history *history, ssize_t j);

#endif /* _HISTORY_H */
