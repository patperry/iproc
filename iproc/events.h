#ifndef _IPROC_EVENTS_H
#define _IPROC_EVENTS_H

#include <stddef.h>
#include <stdint.h>

typedef struct _iproc_events iproc_events;

struct _iproc_events {
    int64_t  ncur;
    int64_t  npast;
    int64_t *cur;
    int64_t *past;
    int64_t *past_dt;
    int64_t  max_ncur;
    int64_t  max_npast;
};

iproc_events *iproc_events_new       ();
void          iproc_events_free      (iproc_events *events);
void          iproc_events_clear     (iproc_events *events);

void          iproc_events_insert    (iproc_events *events,
                                      int64_t       e);
void          iproc_events_advance   (iproc_events *events,
                                      int64_t       dt);

int64_t       iproc_events_ncur      (iproc_events *events);
int64_t       iproc_events_cur       (iproc_events *events,
                                      int64_t       i);
int64_t       iproc_events_find_cur  (iproc_events *events,
                                      int64_t e);

int64_t       iproc_events_npast     (iproc_events *events);
int64_t       iproc_events_past      (iproc_events *events,
                                      int64_t       i);
int64_t       iproc_events_find_past (iproc_events *events,
                                      int64_t e);

int64_t       iproc_events_past_dt   (iproc_events *events,
                                      int64_t       i);


#endif /* _IPROC_EVENTS_H */
