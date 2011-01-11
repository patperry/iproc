#ifndef _IPROC_REFCOUNT_H
#define _IPROC_REFCOUNT_H

#include <iproc/utils.h>

typedef struct _iproc_refcount iproc_refcount;

struct _iproc_refcount {
    int count;
};

void iproc_refcount_init (iproc_refcount *refcount);
void iproc_refcount_set  (iproc_refcount *refcount,
                          int             count);
void iproc_refcount_get  (iproc_refcount *refcount);
int  iproc_refcount_put  (iproc_refcount *refcount,
                          void (*release) (iproc_refcount *refcount));


#endif /* _IPROC_REFCOUNT_H */
