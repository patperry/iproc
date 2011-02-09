#ifndef _IPROC_DESIGN_H
#define _IPROC_DESIGN_H

#include <stdint.h>
#include "actors.h"
#include "history.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"


typedef struct _iproc_design      iproc_design;
typedef struct _iproc_design_ctx  iproc_design_ctx;
typedef struct _iproc_sdesign_var iproc_sdesign_var;


struct _iproc_design {
    iproc_actors  *senders;
    iproc_actors  *receivers;
    int64_t        nstatic;
    int64_t        ndynamic;
    void (*get_sdesign_vars) (iproc_design_ctx *ctx);
    void (*free_user_data) (void *);
    void *user_data;
    iproc_array   *ctxs;
    iproc_refcount refcount;
};

struct _iproc_sdesign_var {
    int64_t        jrecv;
    iproc_svector *jdiff;
};

struct _iproc_design_ctx {
    iproc_design    *design;
    iproc_history *history;
    int64_t        isend;
    iproc_array   *sdesign_vars;
    iproc_refcount refcount;
};

iproc_design *     iproc_design_new           (iproc_actors   *senders,
                                           iproc_actors   *receivers,
                                           int64_t         ndynamic,
                                           void           *user_data,
                                           void          (*get_sdesign_vars) (iproc_design_ctx *ctx),
                                           void          (*free_user_data)  (void *user_data));
iproc_design *     iproc_design_ref           (iproc_design     *design);
void             iproc_design_unref         (iproc_design     *design);

int64_t          iproc_design_dim           (iproc_design     *design);
int64_t          iproc_design_nstatic       (iproc_design     *design);
int64_t          iproc_design_istatic       (iproc_design     *design,
                                           int64_t         i);
int64_t          iproc_design_ndynamic      (iproc_design     *design);
int64_t          iproc_design_idynamic      (iproc_design     *design,
                                           int64_t         i);

int64_t          iproc_design_nsender       (iproc_design     *design);
int64_t          iproc_design_nreceiver     (iproc_design     *design);
iproc_actors *   iproc_design_senders       (iproc_design     *design);
iproc_actors *   iproc_design_receivers     (iproc_design     *design);

void             iproc_design_sender0_mul   (double          alpha,
                                           iproc_trans     trans,
                                           iproc_design     *design,
                                           int64_t         isend,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_vector   *y);
void             iproc_design_sender0_muls  (double          alpha,
                                           iproc_trans     trans,
                                           iproc_design     *design,
                                           int64_t         isend,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_vector   *y);


iproc_design_ctx * iproc_design_ctx_new       (iproc_design     *design,
                                           int64_t         isend,
                                           iproc_history  *h);
iproc_design_ctx * iproc_design_ctx_ref       (iproc_design_ctx *ctx);
void             iproc_design_ctx_unref     (iproc_design_ctx *ctx);


void             iproc_design_ctx_set     (iproc_design_ctx *ctx,
                                           int64_t           isend,
                                           iproc_history    *h);
void             iproc_design_ctx_mul       (double          alpha,
                                           iproc_trans     trans,
                                           iproc_design_ctx *ctx,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_vector   *y);
void             iproc_design_ctx_muls      (double          alpha,
                                           iproc_trans     trans,
                                           iproc_design_ctx *ctx,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_vector   *y);

void             iproc_design_ctx_diff_mul  (double          alpha,
                                           iproc_trans     trans,
                                           iproc_design_ctx *ctx,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_svector  *y);
void             iproc_design_ctx_diff_muls (double          alpha,
                                           iproc_trans     trans,
                                           iproc_design_ctx *ctx,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_svector  *y);


iproc_sdesign_var * iproc_sdesign_var_new  (iproc_design      *design,
                                            int64_t            jrecv);
void                iproc_sdesign_var_free (iproc_design      *design,
                                            iproc_sdesign_var *sv);

#endif /* _IPROC_DESIGN_H */
