#ifndef _IPROC_VARS_H
#define _IPROC_VARS_H

#include <stdint.h>
#include "actors.h"
#include "history.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"


typedef struct _iproc_vars        iproc_vars;
typedef struct _iproc_vars_ctx    iproc_vars_ctx;
typedef struct _iproc_sender_vars iproc_sender_vars;


struct _iproc_vars {
    iproc_actors  *senders;
    iproc_actors  *receivers;
    int64_t        nsender_vars;
    void (*get_sender_vars) (iproc_history *, int64_t, iproc_array *);
    iproc_refcount refcount;
};

struct _iproc_sender_vars {
    int64_t        jrecv;
    iproc_svector *jdiff;
};

struct _iproc_vars_ctx {
    iproc_vars    *vars;
    iproc_history *history;
    int64_t        isend;
    iproc_array   *sender_vars;
    iproc_refcount refcount;
};


iproc_vars *     iproc_vars_new           (iproc_actors   *senders,
                                           iproc_actors   *receivers);
iproc_vars *     iproc_vars_ref           (iproc_vars     *vars);
void             iproc_vars_unref         (iproc_vars     *vars);
int64_t          iproc_vars_dim           (iproc_vars     *vars);
int64_t          iproc_vars_nsender       (iproc_vars     *vars);
int64_t          iproc_vars_nreceiver     (iproc_vars     *vars);
iproc_actors *   iproc_vars_senders       (iproc_vars     *vars);
iproc_actors *   iproc_vars_receivers     (iproc_vars     *vars);

void             iproc_vars_sender0_mul   (double          alpha,
                                           iproc_trans     trans,
                                           iproc_vars     *vars,
                                           int64_t         isend,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_vector   *y);
void             iproc_vars_sender0_muls  (double          alpha,
                                           iproc_trans     trans,
                                           iproc_vars     *vars,
                                           int64_t         isend,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_vector   *y);


iproc_vars_ctx * iproc_vars_ctx_new       (iproc_vars     *vars,
                                           int64_t         isend,
                                           iproc_history  *h);
iproc_vars_ctx * iproc_vars_ctx_ref       (iproc_vars_ctx *ctx);
void             iproc_vars_ctx_unref     (iproc_vars_ctx *ctx);


void             iproc_vars_ctx_mul       (double          alpha,
                                           iproc_trans     trans,
                                           iproc_vars_ctx *ctx,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_vector   *y);
void             iproc_vars_ctx_muls      (double          alpha,
                                           iproc_trans     trans,
                                           iproc_vars_ctx *ctx,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_vector   *y);

void             iproc_vars_ctx_diff_mul  (double          alpha,
                                           iproc_trans     trans,
                                           iproc_vars_ctx *ctx,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_svector  *y);
void             iproc_vars_ctx_diff_muls (double          alpha,
                                           iproc_trans     trans,
                                           iproc_vars_ctx *ctx,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_svector  *y);


#endif /* _IPROC_VARS_H */
