#ifndef _IPROC_PQUEUE_H
#define _IPROC_PQUEUE_H

#include <stdbool.h>
#include <sys/types.h>

#include "darray.h"
#include "compare.h"
#include "refcount.h"

typedef struct _iproc_pqueue iproc_pqueue;

struct _iproc_pqueue {
    struct darray  array;
    compare_fn     compare;
    iproc_refcount refcount;
};

iproc_pqueue * iproc_pqueue_new        (size_t      eltsize,
                                        compare_fn  compare);
iproc_pqueue * iproc_pqueue_new_copy   (iproc_pqueue *pqueue);
void           iproc_pqueue_ref        (iproc_pqueue *pqueue);
void           iproc_pqueue_unref      (iproc_pqueue *pqueue);


bool           iproc_pqueue_empty      (iproc_pqueue *pqueue);
ssize_t        iproc_pqueue_size       (iproc_pqueue *pqueue);
const void *   iproc_pqueue_top        (iproc_pqueue *pqueue);
void *         iproc_pqueue_get_top    (iproc_pqueue *pqueue,
                                        void         *dst);

void *         iproc_pqueue_push       (iproc_pqueue *pqueue,
                                        const void   *eltp);
void *         iproc_pqueue_push_array (iproc_pqueue *pqueue,
                                        const void   *elts,
                                        ssize_t       n);
void           iproc_pqueue_pop        (iproc_pqueue *pqueue);
void *         iproc_pqueue_pop_array  (iproc_pqueue *pqueue,
                                        void         *elts,
                                        ssize_t       n);



#endif /* _IPROC_PQUEUE_H */