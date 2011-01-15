#ifndef _IPROC_SLOGLIK_H
#define _IPROC_SLOGLIK_H

#include "array.h"
#include "history.h"
#include "model.h"
#include "refcount.h"
#include "svector.h"

typedef struct _iproc_sloglik iproc_sloglik;

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

iproc_sloglik * iproc_sloglik_new             (iproc_model   *model,
                                               int64_t        isend);
iproc_sloglik * iproc_sloglik_ref             (iproc_sloglik *sll);
void            iproc_sloglik_unref           (iproc_sloglik *sll);

void            iproc_sloglik_insert          (iproc_sloglik *sll,
                                               iproc_history *history,
                                               int64_t        jrecv);
void            iproc_sloglik_insertm         (iproc_sloglik *sll,
                                               iproc_history *history,
                                               int64_t       *jrecv,
                                               int64_t        n);

double          iproc_vector_acc_sloglik_grad (iproc_vector  *dst_vector,
                                               double         scale,
                                               iproc_sloglik *sll);

#endif /* _IPROC_LOGLIK_H */
