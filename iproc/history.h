#ifndef _IPROC_HISTORY_H
#define _IPROC_HISTORY_H

#include <iproc/events.h>

typedef struct _iproc_history        iproc_history;
typedef struct _iproc_history_events iproc_history_events;

struct _iproc_history_events {
    int64_t       elapsed;
    iproc_events *events;
};

struct _iproc_history {
    int64_t               elapsed;
    iproc_history_events *send;
    iproc_history_events *recv;
    int64_t               nsend;
    int64_t               nrecv;
    int64_t               max_nsend;
    int64_t               max_nrecv;
};


iproc_history * iproc_history_new     ();
void            iproc_history_free    (iproc_history *history);
void            iproc_history_clear   (iproc_history *history);

void            iproc_history_advance (iproc_history *history,
                                       int64_t        dt);
void            iproc_history_insert  (iproc_history *history,
                                       int64_t        from,
                                       int64_t        to);
void            iproc_history_insertm (iproc_history *history,
                                       int            nfrom,
                                       int64_t       *from,
                                       int            nto,
                                       int64_t       *to);

int64_t         iproc_history_nsend   (iproc_history *history);
int64_t         iproc_history_nrecv   (iproc_history *history);
iproc_events *  iproc_history_send    (iproc_history *history,
                                       int64_t        i);
iproc_events *  iproc_history_recv    (iproc_history *history,
                                       int64_t        j);


#endif /* _IPROC_HISTORY_H */
