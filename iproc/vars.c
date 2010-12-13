#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include <iproc/memory.h>
#include <iproc/vars.h>

static void
iproc_vars_free (iproc_vars *vars)
{
    if (vars) {
    }
}

iproc_vars *     
iproc_vars_new (iproc_actors *send,
                iproc_actors *recv)
{
    assert(send);
    assert(recv);

    iproc_vars *vars = iproc_malloc(sizeof(*vars));
    vars->refcount = 1;
    return vars;
}


iproc_vars *
iproc_vars_ref (iproc_vars *vars)
{
    if (vars) {
        vars->refcount = vars->refcount + 1;
    }
    return vars;
}

void
iproc_vars_unref (iproc_vars *vars)
{
    if (!vars)
        return;

    if (vars->refcount == 1) {
        iproc_vars_free(vars);
    } else {
        vars->refcount = vars->refcount - 1;
    }
}

static void
iproc_vars_ctx_free (iproc_vars_ctx *ctx)
{
    if (ctx) {
    }
}

iproc_vars_ctx *
iproc_vars_ctx_new (iproc_vars    *vars,
                    iproc_history *h,
                    int64_t        send)
{
    assert(vars);

    iproc_vars_ctx *ctx = iproc_malloc(sizeof(*ctx));
    ctx->refcount = 1;
    return ctx;
}

iproc_vars_ctx *
iproc_vars_ctx_ref (iproc_vars_ctx *ctx)
{
    if (ctx) {
        ctx->refcount = ctx->refcount + 1;
    }
    return ctx;
}

void
iproc_vars_ctx_unref (iproc_vars_ctx *ctx)
{
    if (!ctx)
        return;

    if (ctx->refcount == 1) {
        iproc_vars_ctx_free(ctx);
    } else {
        ctx->refcount = ctx->refcount - 1;
    }
}

void
iproc_vars_ctx_mul (iproc_trans     trans,
                    double          alpha,
                    iproc_vars_ctx *ctx,
                    iproc_vector   *x,
                    double          beta,
                    iproc_vector   *y)
{
    assert(ctx);
    assert(x);
    assert(y);
}

void
iproc_vars_cxt_diff_mul (iproc_trans     trans,
                         double          alpha,
                         iproc_vars_ctx *ctx,
                         iproc_vector   *x,
                         double          beta,
                         iproc_vector   *y)
{
    assert(ctx);
    assert(x);
    assert(y);
}
