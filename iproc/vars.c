#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include <iproc/vars.h>

iproc_vars *     
iproc_vars_new (iproc_actors *send,
                iproc_actors *recv)
{
    assert(send);
    assert(recv);
    return NULL;
}

void
iproc_vars_free (iproc_vars *vars)
{
    if (vars) {
    }
}

iproc_vars_ctx *
iproc_vars_ctx_new (iproc_vars    *vars,
                    iproc_history *h,
                    int64_t        send)
{
    assert(vars);
    return NULL;
}

void
iproc_vars_ctx_free (iproc_vars_ctx *ctx)
{
    if (ctx) {
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
