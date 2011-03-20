#ifndef _IPROC_TRACE_H
#define _IPROC_TRACE_H

#include <stddef.h>
#include <stdint.h>
#include "array.h"
#include "refcount.h"

typedef struct _iproc_trace iproc_trace;

typedef struct _iproc_event  iproc_event;

struct _iproc_trace {
    double         tcur;
    iproc_array   *cur;
    iproc_array   *past;
    iproc_refcount refcount;
};


struct _iproc_event {
    int64_t event;
    double  t;
};



iproc_trace * iproc_trace_new        ();
iproc_trace * iproc_trace_ref        (iproc_trace *trace);
void          iproc_trace_unref      (iproc_trace *trace);


void          iproc_trace_clear      (iproc_trace *trace);

void          iproc_trace_insert     (iproc_trace *trace,
                                      int64_t       e);
void          iproc_trace_advance_to (iproc_trace *trace,
                                      double        t);

double         iproc_trace_tcur       (iproc_trace *trace);
int64_t        iproc_trace_size       (iproc_trace *trace);
iproc_event * iproc_trace_get        (iproc_trace *trace,
                                        int64_t       i);
iproc_event *  iproc_trace_lookup     (iproc_trace *trace,
                                        int64_t       e);


#endif /* _IPROC_TRACE_H */
