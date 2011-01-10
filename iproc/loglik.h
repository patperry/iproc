#ifndef _IPROC_LOGLIK_H
#define _IPROC_LOGLIK_H

#include <iproc/array.h>
#include <iproc/model.h>
#include <iproc/refcount.h>
#include <iproc/svector.h>

typedef struct _iproc_loglik  iproc_loglik;
typedef struct _iproc_sloglik iproc_sloglik;

struct _iproc_loglik {
    iproc_array *sloglik_array;
};

struct _iproc_sloglik {
    iproc_model   *model;
    int64_t        isend;
    int64_t        nsend;
    int64_t        nrecv;
    iproc_svector *ovarsdiff;
    double         value;
    double         suminvwt;
    double         suminvwt_scale;
    iproc_svector *newprob;
    iproc_svector *evarsdiff;
    iproc_refcount refcount;
};

iproc_loglik * iproc_loglik_new    (iproc_model   *model);
iproc_loglik * iproc_loglik_ref    (iproc_loglik  *loglik);
void           iproc_loglik_unref  (iproc_loglik  *loglik);

void           iproc_loglik_insert (iproc_loglik  *loglik,
                                    int64_t        from,
                                    int64_t        to,
                                    iproc_history *history);

#endif /* _IPROC_LOGLIK_H */
