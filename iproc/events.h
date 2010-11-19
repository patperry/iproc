#ifndef _IPROC_EVENTS_H
#define _IPROC_EVENTS_H

#include <stddef.h>
#include <stdint.h>

typedef struct _iproc_events iproc_events;

struct _iproc_events {
    size_t    ncur;
    size_t    max_ncur;
    size_t    max_npast;
    size_t    npast;
    int64_t  *ecur;
    int64_t  *epast;
    uint64_t *dtpast;
};

iproc_events *iproc_events_new     ();
void          iproc_events_free    (iproc_events *events);
void          iproc_events_clear   (iproc_events *events);

void          iproc_events_insert  (iproc_events *events,
                                    int64_t       e);
void          iproc_events_advance (iproc_events *events,
                                    uint64_t      dt);

size_t        iproc_events_ncur    (iproc_events *events);
int64_t       iproc_events_cur     (iproc_events *events,
                                    size_t        i);

size_t        iproc_events_npast   (iproc_events *events);
int64_t       iproc_events_past    (iproc_events *events,
                                    size_t        i);
uint64_t      iproc_events_dtpast  (iproc_events *events,
                                    size_t        i);

#endif /* _IPROC_EVENTS_H */
