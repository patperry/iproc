#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "design.h"


static void
iproc_design_ctx_free (iproc_design_ctx *ctx)
{
    if (ctx) {
        iproc_history_unref(ctx->history);
        iproc_design_unref(ctx->design);
        
        if (ctx->sender_design) {
            int64_t i, n = iproc_array_size(ctx->sender_design);
            for (i = 0; i < n; i++) {
                iproc_sender_design *sv = &(iproc_array_index(ctx->sender_design,
                                                            iproc_sender_design,
                                                            i));
                iproc_svector_unref(sv->jdiff);
            }
        }

        iproc_array_unref(ctx->sender_design);
        iproc_free(ctx);
    }
}

iproc_design_ctx *
iproc_design_ctx_new (iproc_design    *design,
                    int64_t        isend,
                    iproc_history *h)

{
    assert(design);
    assert(0 <= isend);
    assert(isend < iproc_design_nsender(design));

    iproc_design_ctx *ctx = iproc_malloc(sizeof(*ctx));
    ctx->design = iproc_design_ref(design);
    ctx->history = iproc_history_ref(h);
    ctx->isend = isend;

    iproc_refcount_init(&ctx->refcount);

    if (design->get_sender_design) {
        ctx->sender_design = iproc_array_new(sizeof(iproc_sender_design));
        design->get_sender_design(ctx);
    } else {
        ctx->sender_design = NULL;
    }

    return ctx;
}

iproc_design_ctx *
iproc_design_ctx_ref (iproc_design_ctx *ctx)
{
    if (ctx) {
        iproc_refcount_get(&ctx->refcount);
    }
    return ctx;
}

static void
iproc_design_ctx_release (iproc_refcount *refcount)
{
    iproc_design_ctx *ctx = container_of(refcount, iproc_design_ctx, refcount);
    iproc_design_ctx_free(ctx);
}

void
iproc_design_ctx_unref (iproc_design_ctx *ctx)
{
    if (!ctx)
        return;

    iproc_refcount_put(&ctx->refcount, iproc_design_ctx_release);
}

void
iproc_design_ctx_mul (double          alpha,
                    iproc_trans     trans,
                    iproc_design_ctx *ctx,
                    iproc_vector   *x,
                    double          beta,
                    iproc_vector   *y)
{
    assert(ctx);
    assert(ctx->design);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_design_dim(ctx->design));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_dim(ctx->design));

    iproc_design_sender0_mul(alpha, trans, ctx->design, ctx->isend, x, beta, y);
    
    iproc_svector *diffprod = iproc_svector_new(iproc_vector_dim(y));
    iproc_design_ctx_diff_mul(alpha, trans, ctx, x, 0.0, diffprod);
    iproc_vector_sacc(y, 1.0, diffprod);
    iproc_svector_unref(diffprod);
}

void
iproc_design_ctx_muls (double          alpha,
                     iproc_trans     trans,
                     iproc_design_ctx *ctx,
                     iproc_svector  *x,
                     double          beta,
                     iproc_vector   *y)
{
    assert(ctx);
    assert(ctx->design);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_design_dim(ctx->design));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_design_dim(ctx->design));

    iproc_design_sender0_muls(alpha, trans, ctx->design, ctx->isend, x, beta, y);
    
    iproc_svector *diffprod = iproc_svector_new(iproc_vector_dim(y));
    iproc_design_ctx_diff_muls(alpha, trans, ctx, x, 0.0, diffprod);
    iproc_vector_sacc(y, 1.0, diffprod);
    iproc_svector_unref(diffprod);
}


void
iproc_design_ctx_diff_mul (double          alpha,
                         iproc_trans     trans,
                         iproc_design_ctx *ctx,
                         iproc_vector   *x,
                         double          beta,
                         iproc_svector  *y)
{
    assert(ctx);
    assert(ctx->design);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_design_dim(ctx->design));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_design_dim(ctx->design));

    
    /* y := beta y */
    if (beta == 0.0) {
        iproc_svector_clear(y);
    } else if (beta != 1.0) {
        iproc_svector_scale(y, beta);
    }

    if (ctx->history == NULL)
        return;

    iproc_design *design = ctx->design;
    int64_t ndynamic = iproc_design_ndynamic(design);

    if (ndynamic == 0)
        return;
    
    int64_t ix_begin = iproc_design_idynamic(design, 0);
    int64_t ix_end = ix_begin + ndynamic;
    iproc_array *sender_design = ctx->sender_design;

    if (trans == IPROC_TRANS_NOTRANS) {
        iproc_vector_view xsub = iproc_vector_subvector(x, ix_begin, ndynamic);
        int64_t i, n = iproc_array_size(sender_design);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_design *sv = &(iproc_array_index(sender_design,
                                                        iproc_sender_design,
                                                        i));
            int64_t jrecv = sv->jrecv;
            iproc_svector *jdiff = sv->jdiff;
            double dot = iproc_vector_sdot(&xsub.vector, jdiff);
            iproc_svector_inc(y, jrecv, alpha * dot);
        }
    } else {
        int64_t i, n = iproc_array_size(sender_design);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_design *sv = &(iproc_array_index(sender_design,
                                                        iproc_sender_design,
                                                        i));
            int64_t jrecv = sv->jrecv;
            iproc_svector *jdiff = sv->jdiff;
            double xjrecv = iproc_vector_get(x, jrecv);

            if (xjrecv == 0.0)
                continue;

            /* y := y + alpha * x[j] * diff[j] */
            double jscale = alpha * xjrecv;
            int64_t inz, nnz = iproc_svector_nnz(jdiff);
            for (inz = 0; inz < nnz; inz++) {
                int64_t ix = iproc_svector_nz(jdiff, inz);
                double val = iproc_svector_nz_get(jdiff, inz);

                assert(ix + ix_begin < ix_end);
                iproc_svector_inc(y, ix + ix_begin, jscale * val);
            }
        }
    }
}


void
iproc_design_ctx_diff_muls (double          alpha,
                          iproc_trans     trans,
                          iproc_design_ctx *ctx,
                          iproc_svector  *x,
                          double          beta,
                          iproc_svector  *y)
{
    assert(ctx);
    assert(ctx->design);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_design_dim(ctx->design));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_design_nreceiver(ctx->design));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_design_dim(ctx->design));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_svector_clear(y);
    } else if (beta != 1.0) {
        iproc_svector_scale(y, beta);
    }

    if (ctx->history == NULL)
        return;

    iproc_design *design = ctx->design;
    int64_t ndynamic = iproc_design_ndynamic(design);

    if (ndynamic == 0)
        return;
    
    int64_t ix_begin = iproc_design_idynamic(design, 0);
    int64_t ix_end = ix_begin + ndynamic;
    iproc_array *sender_design = ctx->sender_design;

    if (trans == IPROC_TRANS_NOTRANS) {
        int64_t i, n = iproc_array_size(sender_design);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_design *sv = &(iproc_array_index(sender_design,
                                                        iproc_sender_design,
                                                        i));
            int64_t jrecv = sv->jrecv;
            iproc_svector *jdiff = sv->jdiff;
            int64_t inz, nnz = iproc_svector_nnz(x);
            double dot = 0.0; 
            for (inz = 0; inz < nnz; inz++) {
                int64_t ix = iproc_svector_nz(x, inz);
                if (ix < ix_begin || ix >= ix_end)
                    continue;

                double xval = iproc_svector_nz_get(x, inz);
                double diffval = iproc_svector_get(jdiff, ix - ix_begin);
                dot += xval * diffval;
            }

            iproc_svector_inc(y, jrecv, alpha * dot);
        }
    } else {
        int64_t i, n = iproc_array_size(sender_design);
        
        for (i = 0; i < n; i++) {
            iproc_sender_design *sv = &(iproc_array_index(sender_design,
                                                        iproc_sender_design,
                                                        i));
            int64_t jrecv = sv->jrecv;
            iproc_svector *jdiff = sv->jdiff;
            double xjrecv = iproc_svector_get(x, jrecv);
            if (xjrecv == 0.0)
                continue;

            /* y := y + alpha * x[j] * diff[j] */
            double jscale = alpha * xjrecv;
            int64_t inz, nnz = iproc_svector_nnz(jdiff);
            for (inz = 0; inz < nnz; inz++) {
                int64_t ix = iproc_svector_nz(jdiff, inz);
                double val = iproc_svector_nz_get(jdiff, inz);

                assert(ix + ix_begin < ix_end);
                iproc_svector_inc(y, ix + ix_begin, jscale * val);
            }
        }
    }

}
