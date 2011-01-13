#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include "memory.h"
#include "vars.h"

static void
iproc_vars_free (iproc_vars *vars)
{
    if (vars) {
        iproc_actors_unref(vars->receivers);
        iproc_actors_unref(vars->senders);
        iproc_free(vars);
    }
}

iproc_vars *     
iproc_vars_new (iproc_actors *senders,
                iproc_actors *receivers)
{
    assert(senders);
    assert(receivers);

    iproc_vars *vars = iproc_malloc(sizeof(*vars));
    vars->senders = iproc_actors_ref(senders);
    vars->receivers = iproc_actors_ref(receivers);
    iproc_refcount_init(&vars->refcount);
    return vars;
}

iproc_vars *
iproc_vars_ref (iproc_vars *vars)
{
    if (vars) {
        iproc_refcount_get(&vars->refcount);
    }
    return vars;
}

static void
iproc_vars_release (iproc_refcount *refcount)
{
    iproc_vars *vars = container_of(refcount, iproc_vars, refcount);
    iproc_vars_free(vars);
}

void
iproc_vars_unref (iproc_vars *vars)
{
    if (!vars)
        return;

    iproc_refcount_put(&vars->refcount, iproc_vars_release);
}

static void
iproc_vars_ctx_free (iproc_vars_ctx *ctx)
{
    if (ctx) {
        iproc_history_unref(ctx->history);
        iproc_vars_unref(ctx->vars);
        iproc_free(ctx);
    }
}

iproc_vars_ctx *
iproc_vars_ctx_new (iproc_vars    *vars,
                    iproc_history *h,
                    int64_t        isend)
{
    assert(vars);
    assert(0 <= isend);
    assert(isend < iproc_vars_nsender(vars));

    iproc_vars_ctx *ctx = iproc_malloc(sizeof(*ctx));
    ctx->vars = iproc_vars_ref(vars);
    ctx->history = iproc_history_ref(h);
    ctx->isend = isend;
    iproc_refcount_init(&ctx->refcount);
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

int64_t
iproc_vars_dim (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *senders = iproc_vars_senders(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    int64_t dim = p * q;
    return dim;
}

int64_t
iproc_vars_nsender (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *senders = iproc_vars_senders(vars);
    return iproc_actors_size(senders);
}

int64_t
iproc_vars_nreceiver (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    return iproc_actors_size(receivers);
}

iproc_actors *
iproc_vars_senders (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *senders = vars->senders;
    return senders;
}

iproc_actors *
iproc_vars_receivers (iproc_vars *vars)
{
    assert(vars);
    iproc_actors *receivers = vars->receivers;
    return receivers;
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

    iproc_vars *vars = ctx->vars;
    iproc_actors *senders = iproc_vars_senders(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    int64_t p = iproc_actors_dim(senders);
    int64_t q = iproc_actors_dim(receivers);
    iproc_vector *s = iproc_actors_traits(senders, ctx->isend);
    iproc_vector *z = iproc_vector_new(q);

    /* y := beta y */
    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (trans == IPROC_TRANS_NOTRANS) {
        iproc_vector_view xsub = iproc_vector_subvector(x, 0, p * q);

        /* z := alpha t(x) s */
        iproc_matrix_view xmat = iproc_matrix_view_vector(&xsub.vector, p, q);
        iproc_matrix_mul(alpha, IPROC_TRANS_TRANS, &xmat.matrix, s, 0.0, z);

        /* y := y + R z */
        iproc_actors_mul(1.0, IPROC_TRANS_NOTRANS, receivers, z, 1.0, y);
    } else {
        /* z := alpha t(R) x */
        iproc_actors_mul(alpha, IPROC_TRANS_TRANS, receivers, x, 0.0, z);

        /* y := y + s \otimes z */
        iproc_vector_view ysub = iproc_vector_subvector(y, 0, p * q);
        iproc_matrix_view ymat = iproc_matrix_view_vector(&ysub.vector, p, q);
        iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
        iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
        iproc_matrix_matmul(1.0, IPROC_TRANS_NOTRANS, &smat.matrix, &zmat.matrix,
                            1.0, &ymat.matrix);
    }

    /* TODO: need something like the following:
     * iproc_vars_ctx_diff_mul(alpha, trans, ctx, x, 1.0, y);
     */
    iproc_vector_unref(z);
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
}
