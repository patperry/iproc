#ifndef _IPROC_EVENTS_H
#define _IPROC_EVENTS_H

#include <stddef.h>
#include <stdint.h>
#include <iproc/array.h>

typedef struct _iproc_event      iproc_event;
typedef struct _iproc_events     iproc_events;
typedef struct _iproc_past_event iproc_past_event;

struct _iproc_events {
    iproc_array *cur;
    iproc_array *past;
};

struct _iproc_event {
    int64_t e;
};

struct _iproc_past_event {
    iproc_event event;
    int64_t     dt;
};


iproc_events * iproc_events_new       ();
void           iproc_events_free      (iproc_events *events);
void           iproc_events_clear     (iproc_events *events);

void           iproc_events_insert    (iproc_events *events,
                                       int64_t       e);
void           iproc_events_advance   (iproc_events *events,
                                       int64_t       dt);

int64_t        iproc_events_ncur      (iproc_events *events);
int64_t        iproc_events_cur       (iproc_events *events,
                                       int64_t       i);
int64_t        iproc_events_find_cur  (iproc_events *events,
                                       int64_t       e);

int64_t        iproc_events_npast     (iproc_events *events);
int64_t        iproc_events_past      (iproc_events *events,
                                       int64_t       i);
int64_t        iproc_events_find_past (iproc_events *events,
                                       int64_t       e);

int64_t        iproc_events_past_dt   (iproc_events *events,
                                       int64_t       i);


#endif /* _IPROC_EVENTS_H */
