#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "vars.h"


static void
iproc_vars_ctx_free (iproc_vars_ctx *ctx)
{
    if (ctx) {
        iproc_history_unref(ctx->history);
        iproc_vars_unref(ctx->vars);
        
        if (ctx->sender_vars) {
            int64_t i, n = iproc_array_size(ctx->sender_vars);
            for (i = 0; i < n; i++) {
                iproc_sender_vars *sv = &(iproc_array_index(ctx->sender_vars,
                                                            iproc_sender_vars,
                                                            i));
                iproc_svector_unref(sv->jdiff);
            }
        }

        iproc_array_unref(ctx->sender_vars);
        iproc_free(ctx);
    }
}

iproc_vars_ctx *
iproc_vars_ctx_new (iproc_vars    *vars,
                    int64_t        isend,
                    iproc_history *h)

{
    assert(vars);
    assert(0 <= isend);
    assert(isend < iproc_vars_nsender(vars));

    iproc_vars_ctx *ctx = iproc_malloc(sizeof(*ctx));
    ctx->vars = iproc_vars_ref(vars);
    ctx->history = iproc_history_ref(h);
    ctx->isend = isend;

    iproc_refcount_init(&ctx->refcount);

    if (vars->get_sender_vars) {
        ctx->sender_vars = iproc_array_new(sizeof(iproc_sender_vars));
        vars->get_sender_vars(ctx);
    } else {
        ctx->sender_vars = NULL;
    }

    return ctx;
}

iproc_vars_ctx *
iproc_vars_ctx_ref (iproc_vars_ctx *ctx)
{
    if (ctx) {
        iproc_refcount_get(&ctx->refcount);
    }
    return ctx;
}

static void
iproc_vars_ctx_release (iproc_refcount *refcount)
{
    iproc_vars_ctx *ctx = container_of(refcount, iproc_vars_ctx, refcount);
    iproc_vars_ctx_free(ctx);
}

void
iproc_vars_ctx_unref (iproc_vars_ctx *ctx)
{
    if (!ctx)
        return;

    iproc_refcount_put(&ctx->refcount, iproc_vars_ctx_release);
}

void
iproc_vars_ctx_mul (double          alpha,
                    iproc_trans     trans,
                    iproc_vars_ctx *ctx,
                    iproc_vector   *x,
                    double          beta,
                    iproc_vector   *y)
{
    assert(ctx);
    assert(ctx->vars);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_dim(ctx->vars));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_dim(ctx->vars));

    iproc_vars_sender0_mul(alpha, trans, ctx->vars, ctx->isend, x, beta, y);
    
    iproc_svector *diffprod = iproc_svector_new(iproc_vector_dim(y));
    iproc_vars_ctx_diff_mul(alpha, trans, ctx, x, 0.0, diffprod);
    iproc_vector_sacc(y, 1.0, diffprod);
    iproc_svector_unref(diffprod);
}

void
iproc_vars_ctx_muls (double          alpha,
                     iproc_trans     trans,
                     iproc_vars_ctx *ctx,
                     iproc_svector  *x,
                     double          beta,
                     iproc_vector   *y)
{
    assert(ctx);
    assert(ctx->vars);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_vars_dim(ctx->vars));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_dim(ctx->vars));

    iproc_vars_sender0_muls(alpha, trans, ctx->vars, ctx->isend, x, beta, y);
    
    iproc_svector *diffprod = iproc_svector_new(iproc_vector_dim(y));
    iproc_vars_ctx_diff_muls(alpha, trans, ctx, x, 0.0, diffprod);
    iproc_vector_sacc(y, 1.0, diffprod);
    iproc_svector_unref(diffprod);
}


void
iproc_vars_ctx_diff_mul (double          alpha,
                         iproc_trans     trans,
                         iproc_vars_ctx *ctx,
                         iproc_vector   *x,
                         double          beta,
                         iproc_svector  *y)
{
    assert(ctx);
    assert(ctx->vars);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_dim(ctx->vars));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_vars_dim(ctx->vars));

    
    /* y := beta y */
    if (beta == 0.0) {
        iproc_svector_clear(y);
    } else if (beta != 1.0) {
        iproc_svector_scale(y, beta);
    }

    if (ctx->history == NULL)
        return;

    iproc_vars *vars = ctx->vars;
    int64_t ndynamic = iproc_vars_ndynamic(vars);

    if (ndynamic == 0)
        return;
    
    int64_t ix_begin = iproc_vars_idynamic(vars, 0);
    int64_t ix_end = ix_begin + ndynamic;
    iproc_array *sender_vars = ctx->sender_vars;

    if (trans == IPROC_TRANS_NOTRANS) {
        iproc_vector_view xsub = iproc_vector_subvector(x, ix_begin, ndynamic);
        int64_t i, n = iproc_array_size(sender_vars);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_vars *sv = &(iproc_array_index(sender_vars,
                                                        iproc_sender_vars,
                                                        i));
            int64_t jrecv = sv->jrecv;
            iproc_svector *jdiff = sv->jdiff;
            double dot = iproc_vector_sdot(&xsub.vector, jdiff);
            iproc_svector_inc(y, jrecv, alpha * dot);
        }
    } else {
        int64_t i, n = iproc_array_size(sender_vars);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_vars *sv = &(iproc_array_index(sender_vars,
                                                        iproc_sender_vars,
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
                double val = iproc_svector_nz_val(jdiff, inz);

                assert(ix + ix_begin < ix_end);
                iproc_svector_inc(y, ix + ix_begin, jscale * val);
            }
        }
    }
}


void
iproc_vars_ctx_diff_muls (double          alpha,
                          iproc_trans     trans,
                          iproc_vars_ctx *ctx,
                          iproc_svector  *x,
                          double          beta,
                          iproc_svector  *y)
{
    assert(ctx);
    assert(ctx->vars);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_vars_dim(ctx->vars));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_vars_nreceiver(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(y) == iproc_vars_dim(ctx->vars));

    /* y := beta y */
    if (beta == 0.0) {
        iproc_svector_clear(y);
    } else if (beta != 1.0) {
        iproc_svector_scale(y, beta);
    }

    if (ctx->history == NULL)
        return;

    iproc_vars *vars = ctx->vars;
    int64_t ndynamic = iproc_vars_ndynamic(vars);

    if (ndynamic == 0)
        return;
    
    int64_t ix_begin = iproc_vars_idynamic(vars, 0);
    int64_t ix_end = ix_begin + ndynamic;
    iproc_array *sender_vars = ctx->sender_vars;

    if (trans == IPROC_TRANS_NOTRANS) {
        int64_t i, n = iproc_array_size(sender_vars);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_vars *sv = &(iproc_array_index(sender_vars,
                                                        iproc_sender_vars,
                                                        i));
            int64_t jrecv = sv->jrecv;
            iproc_svector *jdiff = sv->jdiff;
            int64_t inz, nnz = iproc_svector_nnz(x);
            double dot = 0.0; 
            for (inz = 0; inz < nnz; inz++) {
                int64_t ix = iproc_svector_nz(x, inz);
                if (ix < ix_begin || ix >= ix_end)
                    continue;

                double xval = iproc_svector_nz_val(x, inz);
                double diffval = iproc_svector_get(jdiff, ix - ix_begin);
                dot += xval * diffval;
            }

            iproc_svector_inc(y, jrecv, alpha * dot);
        }
    } else {
        int64_t i, n = iproc_array_size(sender_vars);
        
        for (i = 0; i < n; i++ ) {
            iproc_sender_vars *sv = &(iproc_array_index(sender_vars,
                                                        iproc_sender_vars,
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
                double val = iproc_svector_nz_val(jdiff, inz);

                assert(ix + ix_begin < ix_end);
                iproc_svector_inc(y, ix + ix_begin, jscale * val);
            }
        }
    }

}
