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

