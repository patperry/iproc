#ifndef _IPROC_HISTORY_H
#define _IPROC_HISTORY_H

#include "array.h"
#include "events.h"
#include "refcount.h"

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
