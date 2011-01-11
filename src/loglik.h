#ifndef _IPROC_LOGLIK_H
#define _IPROC_LOGLIK_H

#include "array.h"
#include "model.h"
#include "refcount.h"
#include "svector.h"

typedef struct _iproc_loglik iproc_loglik;

struct _iproc_loglik {
    iproc_model   *model;
    iproc_array   *sloglik_array;
    iproc_refcount refcount;
};


iproc_loglik * iproc_loglik_new             (iproc_model   *model);
iproc_loglik * iproc_loglik_ref             (iproc_loglik  *loglik);
void           iproc_loglik_unref           (iproc_loglik  *loglik);

void           iproc_loglik_insert          (iproc_loglik  *loglik,
                                             int64_t        from,
                                             int64_t        to,
                                             iproc_history *history);

double         iproc_vector_acc_loglik_grad (iproc_vector *dst_vector,
                                             double        scale,
                                             iproc_loglik *loglik);

#endif /* _IPROC_LOGLIK_H */
