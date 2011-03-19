#ifndef _IPROC_HISTORY_H
#define _IPROC_HISTORY_H

#include "array.h"
#include "events.h"
#include "refcount.h"


/* A `history` object stores the history of the process up to the
 * current time.  Internally, the history maintains two arrays, one
 * containing the sender events, and one containing the receiver events.
 * The message (i -> j) gets recorded as two events:
 *
 *         (1) event "j" in send[i]; and
 *         (2) event "i" in recv[j].
 *
 * Events can be added to the history by `iproc_history_insert` and
 * `iproc_history_insertm`, the latter of which is for multicast events.
 * The current time can be changed by calling `iiproc_history_advance_to`.
 *
 * The function `iproc_history_send(h,i)` gets the events `send[i]`; the
 * function `iiproc_hiostyr_recv(h,j)` gets the events `recv[j]`.
 */

typedef struct _iproc_history        iproc_history;
typedef struct _iproc_history_events iproc_history_events;

struct _iproc_history_events {
    iproc_events *events;
    double        tcur;
};

struct _iproc_history {
    double         tcur;
    iproc_array   *send;
    iproc_array   *recv;
    iproc_refcount refcount;
};


iproc_history * iproc_history_new     ();
iproc_history * iproc_history_ref     (iproc_history *history);
void            iproc_history_unref   (iproc_history *history);
void            iproc_history_clear   (iproc_history *history);

double          iproc_hisotry_tcur    (iproc_history *history);
void            iproc_history_advance_to (iproc_history *history,
                                          double         t);
void            iproc_history_insert  (iproc_history *history,
                                       int64_t        from,
                                       int64_t        to);
void            iproc_history_insertm (iproc_history *history,
                                       int64_t        from,
                                       int64_t       *to,
                                       int64_t        nto);

int64_t         iproc_history_nsend   (iproc_history *history);
int64_t         iproc_history_nrecv   (iproc_history *history);
iproc_events *  iproc_history_send    (iproc_history *history,
                                       int64_t        i);
iproc_events *  iproc_history_recv    (iproc_history *history,
                                       int64_t        j);


#endif /* _IPROC_HISTORY_H */
