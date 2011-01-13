#ifndef _IPROC_VARS_H
#define _IPROC_VARS_H

#include <stdint.h>
#include "actors.h"
#include "history.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"

typedef struct _iproc_vars     iproc_vars;
typedef struct _iproc_vars_ctx iproc_vars_ctx;


struct _iproc_vars {
    iproc_actors  *senders;
    iproc_actors  *receivers;
    iproc_refcount refcount;
};

struct _iproc_vars_ctx {
    iproc_vars    *vars;
    iproc_history *history;
    int64_t        isend;
    iproc_refcount refcount;
};


iproc_vars *     iproc_vars_new          (iproc_actors   *senders,
                                          iproc_actors   *receivers);
iproc_vars *     iproc_vars_ref          (iproc_vars     *vars);
void             iproc_vars_unref        (iproc_vars     *vars);
int64_t          iproc_vars_dim          (iproc_vars     *vars);
int64_t          iproc_vars_nsender      (iproc_vars     *vars);
int64_t          iproc_vars_nreceiver    (iproc_vars     *vars);
iproc_actors *   iproc_vars_senders      (iproc_vars     *vars);
iproc_actors *   iproc_vars_receivers    (iproc_vars     *vars);

iproc_vars_ctx * iproc_vars_ctx_new      (iproc_vars     *vars,
                                          iproc_history  *h,
                                          int64_t         isend);
iproc_vars_ctx * iproc_vars_ctx_ref      (iproc_vars_ctx *ctx);
void             iproc_vars_ctx_unref    (iproc_vars_ctx *ctx);

void             iproc_vars_ctx_mul      (double          alpha,
                                          iproc_trans     trans,
                                          iproc_vars_ctx *ctx,
                                          iproc_vector   *x,
                                          double          beta,
                                          iproc_vector   *y);
void             iproc_vars_ctx_diff_mul (double          alpha,
                                          iproc_trans     trans,
                                          iproc_vars_ctx *ctx,
                                          iproc_vector   *x,
                                          double          beta,
                                          iproc_svector  *y);

#endif /* _IPROC_VARS_H */
