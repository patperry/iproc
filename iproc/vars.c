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
        iproc_actors_unref(vars->recv);
        iproc_actors_unref(vars->send);
        iproc_free(vars);
    }
}

iproc_vars *     
iproc_vars_new (iproc_actors *send,
                iproc_actors *recv)
{
    assert(send);
    assert(recv);

    iproc_vars *vars = iproc_malloc(sizeof(*vars));
    vars->send = iproc_actors_ref(send);
    vars->recv = iproc_actors_ref(recv);
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
    assert(vars->send);
    assert(0 <= isend);
    assert(isend < iproc_actors_size(vars->send));

    iproc_vars_ctx *ctx = iproc_malloc(sizeof(*ctx));
    ctx->vars = iproc_vars_ref(vars);
    ctx->history = iproc_history_ref(h);
    ctx->isend = isend;
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

int64_t
iproc_vars_dim (iproc_vars *vars)
{
    assert(vars);
    int64_t p = iproc_actors_dim(vars->send);
    int64_t q = iproc_actors_dim(vars->recv);
    int64_t dim = p * q;
    return dim;
}

int64_t
iproc_vars_nsend (iproc_vars *vars)
{
    assert(vars);
    return iproc_actors_size(vars->send);
}

int64_t
iproc_vars_nrecv (iproc_vars *vars)
{
    assert(vars);
    return iproc_actors_size(vars->recv);
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
           || iproc_vector_dim(y) == iproc_vars_nrecv(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_nrecv(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_dim(ctx->vars));

    iproc_vars *vars = ctx->vars;
    int64_t p = iproc_actors_dim(vars->send);
    int64_t q = iproc_actors_dim(vars->recv);
    iproc_vector *s = iproc_actors_vector(vars->send, ctx->isend);
    iproc_vector *z = iproc_vector_new(q);

    if (trans == IPROC_TRANS_NOTRANS) {
        /* z := alpha t(x) s */
        iproc_matrix_view xmat = iproc_matrix_view_vector(x, p, q);
        iproc_matrix_mul(alpha, IPROC_TRANS_TRANS, &xmat.matrix, s, 0.0, z);

        /* y := beta y + R z */
        iproc_actors_mul(1.0, IPROC_TRANS_NOTRANS, vars->recv, z, beta, y);
    } else {
        /* z := alpha t(R) x */
        iproc_actors_mul(alpha, IPROC_TRANS_TRANS, vars->recv, x, 0.0, z);

        /* y := beta y + s \otimes z */
        iproc_matrix_view ymat = iproc_matrix_view_vector(y, p, q);
        iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
        iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
        iproc_matrix_matmul(1.0, IPROC_TRANS_NOTRANS, &smat.matrix, &zmat.matrix,
                            beta, &ymat.matrix);
    }

    iproc_vars_ctx_diff_mul(alpha, trans, ctx, x, 1.0, y);
    iproc_vector_unref(z);
}

void
iproc_vars_ctx_diff_mul (double          alpha,
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
           || iproc_vector_dim(y) == iproc_vars_nrecv(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_vars_nrecv(ctx->vars));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_vars_dim(ctx->vars));

    if (beta == 0.0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1.0) {
        iproc_vector_scale(y, beta);
    }

    if (ctx->history == NULL)
        return;
}
