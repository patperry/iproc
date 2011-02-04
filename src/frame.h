#ifndef _IPROC_FRAME_H
#define _IPROC_FRAME_H

#include <stdint.h>
#include "actors.h"
#include "history.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"


typedef struct _iproc_frame        iproc_frame;
typedef struct _iproc_frame_ctx    iproc_frame_ctx;
typedef struct _iproc_sender_frame iproc_sender_frame;


struct _iproc_frame {
    iproc_actors  *senders;
    iproc_actors  *receivers;
    int64_t        nstatic;
    int64_t        ndynamic;
    void (*get_sender_frame) (iproc_frame_ctx *ctx);
    void (*free_user_data) (void *);
    void *user_data;
    iproc_refcount refcount;
};

struct _iproc_sender_frame {
    int64_t        jrecv;
    iproc_svector *jdiff;
};

struct _iproc_frame_ctx {
    iproc_frame    *frame;
    iproc_history *history;
    int64_t        isend;
    iproc_array   *sender_frame;
    iproc_refcount refcount;
};

iproc_frame *     iproc_frame_new           (iproc_actors   *senders,
                                           iproc_actors   *receivers,
                                           int64_t         ndynamic,
                                           void           *user_data,
                                           void          (*get_sender_frame) (iproc_frame_ctx *ctx),
                                           void          (*free_user_data)  (void *user_data));
iproc_frame *     iproc_frame_ref           (iproc_frame     *frame);
void             iproc_frame_unref         (iproc_frame     *frame);

int64_t          iproc_frame_dim           (iproc_frame     *frame);
int64_t          iproc_frame_nstatic       (iproc_frame     *frame);
int64_t          iproc_frame_istatic       (iproc_frame     *frame,
                                           int64_t         i);
int64_t          iproc_frame_ndynamic      (iproc_frame     *frame);
int64_t          iproc_frame_idynamic      (iproc_frame     *frame,
                                           int64_t         i);

int64_t          iproc_frame_nsender       (iproc_frame     *frame);
int64_t          iproc_frame_nreceiver     (iproc_frame     *frame);
iproc_actors *   iproc_frame_senders       (iproc_frame     *frame);
iproc_actors *   iproc_frame_receivers     (iproc_frame     *frame);

void             iproc_frame_sender0_mul   (double          alpha,
                                           iproc_trans     trans,
                                           iproc_frame     *frame,
                                           int64_t         isend,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_vector   *y);
void             iproc_frame_sender0_muls  (double          alpha,
                                           iproc_trans     trans,
                                           iproc_frame     *frame,
                                           int64_t         isend,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_vector   *y);


iproc_frame_ctx * iproc_frame_ctx_new       (iproc_frame     *frame,
                                           int64_t         isend,
                                           iproc_history  *h);
iproc_frame_ctx * iproc_frame_ctx_ref       (iproc_frame_ctx *ctx);
void             iproc_frame_ctx_unref     (iproc_frame_ctx *ctx);


void             iproc_frame_ctx_mul       (double          alpha,
                                           iproc_trans     trans,
                                           iproc_frame_ctx *ctx,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_vector   *y);
void             iproc_frame_ctx_muls      (double          alpha,
                                           iproc_trans     trans,
                                           iproc_frame_ctx *ctx,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_vector   *y);

void             iproc_frame_ctx_diff_mul  (double          alpha,
                                           iproc_trans     trans,
                                           iproc_frame_ctx *ctx,
                                           iproc_vector   *x,
                                           double          beta,
                                           iproc_svector  *y);
void             iproc_frame_ctx_diff_muls (double          alpha,
                                           iproc_trans     trans,
                                           iproc_frame_ctx *ctx,
                                           iproc_svector  *x,
                                           double          beta,
                                           iproc_svector  *y);


#endif /* _IPROC_FRAME_H */
