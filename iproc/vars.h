#ifndef _IPROC_VARS_H
#define _IPROC_VARS_H

#include <stdint.h>
#include <iproc/actors.h>
#include <iproc/blas.h>
#include <iproc/history.h>
#include <iproc/vector.h>

typedef struct _iproc_vars     iproc_vars;
typedef struct _iproc_vars_ctx iproc_vars_ctx;


struct _iproc_vars {
    iproc_actors *send;
    iproc_actors *recv;
    int           refcount;
};

struct _iproc_vars_ctx {
    iproc_vars    *vars;
    iproc_history *history;
    int64_t        send;
    int            refcount;
};


iproc_vars *     iproc_vars_new          (iproc_actors   *send,
                                          iproc_actors   *recv);
iproc_vars *     iproc_vars_ref          (iproc_vars     *vars);
void             iproc_vars_unref        (iproc_vars     *vars);
iproc_vars_ctx * iproc_vars_ctx_new      (iproc_vars     *vars,
                                          iproc_history  *h,
                                          int64_t         send);
iproc_vars_ctx * iproc_vars_ctx_ref      (iproc_vars_ctx *ctx);
void             iproc_vars_ctx_unref    (iproc_vars_ctx *ctx);

void             iproc_vars_ctx_mul      (iproc_trans     trans,
                                          double          alpha,
                                          iproc_vars_ctx *ctx,
                                          iproc_vector   *x,
                                          double          beta,
                                          iproc_vector   *y);

void             iproc_vars_cxt_diff_mul (iproc_trans     trans,
                                          double          alpha,
                                          iproc_vars_ctx *ctx,
                                          iproc_vector   *x,
                                          double          beta,
                                          iproc_vector   *y);

#endif /* _IPROC_VARS_H */
