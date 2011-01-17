#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "logsumexp.h"
#include "memory.h"
#include "refcount.h"
#include "model.h"

/* The output of this function (logprob0_shift and logprobs) can be used
 * to compute the probabilites for the given context.  To compute log(p[j])
 * do the following:
 *
 *   1. If jrecv is a sparse index of logprobs, return logprobs[jrecv]
 *   2. Otherwise return logprobs0[jrecv] - logsumweight_diff
 *
 */
static void
compute_new_logprobs (iproc_vars_ctx *ctx,
                      iproc_vector   *coefs,
                      int             has_loops,
                      iproc_vector   *logprobs0,
                      iproc_svector  *logprobs,
                      double         *plogsumweight_diff)
{
    assert(ctx);
    assert(coefs);
    assert(logprobs0);
    assert(logprobs);
    assert(plogsumweight_diff);
    assert(iproc_vector_dim(logprobs0) == iproc_svector_dim(logprobs));

    /* compute the changes in weights */
    iproc_vars_ctx_diff_mul(1.0, IPROC_TRANS_NOTRANS, ctx, coefs, 0.0, logprobs);
    if (!has_loops) {
        int64_t isend = ctx->isend;
        iproc_svector_set(logprobs, isend, -INFINITY);
    }

    int64_t i, nnz = iproc_svector_nnz(logprobs);
    iproc_vector_view logprobs_nz = iproc_svector_view_nz(logprobs);

    /* compute the scale for the weight differences */
    double lwmax = iproc_vector_max(&logprobs_nz.vector);
    double logscale = IPROC_MAX(0.0, lwmax);
    double invscale = exp(-logscale);

    iproc_logsumexp pos, neg;
    iproc_logsumexp_init(&pos);
    iproc_logsumexp_init(&neg);

    /* treat the initial sum of weights as the first positive difference */
    iproc_logsumexp_insert(&pos, -logscale);

    /* compute the log sums of the positive and negative differences in weights */
    for (i = 0; i < nnz; i++) {
        int64_t jrecv = iproc_svector_nz(logprobs, i);
        double lp0 = iproc_vector_get(logprobs0, jrecv);
        double dlw = iproc_vector_get(&logprobs_nz.vector, i);
        double log_abs_dw;

        /* When w > w0:
         *   w - w0      = w0 * [exp(dlw) - 1];
         *               = w0 * {exp[dlw - log(scale)] - 1/scale} * scale
         *
         *   log(w - w0) = log(w0) + log{exp[(dlw) - log(scale)] - 1/scale} + log(scale)
         */
        if (dlw >= 0) {
            log_abs_dw = lp0 + log(exp(dlw - logscale) - invscale);
            iproc_logsumexp_insert(&pos, log_abs_dw);
        } else {
            log_abs_dw = lp0 + log(invscale - exp(dlw - logscale));
            iproc_logsumexp_insert(&neg, log_abs_dw);
        }

        /* make logprobs[j] store the (unnormalized) log probability */
        iproc_vector_inc(&logprobs_nz.vector, i, lp0);
    }

    double log_sum_abs_dw_p = iproc_logsumexp_value(&pos);
    double log_sum_abs_dw_n = iproc_logsumexp_value(&neg);

    if (log_sum_abs_dw_n > log_sum_abs_dw_p) {
        /* The sum of the weights is positive, so this only happens as a
         * result of numerical errors */
        log_sum_abs_dw_n = log_sum_abs_dw_p;
        printf("\nWARNING: numerical errors in model-ctx.c (compute_new_logprobs)"
               "\nPlease report a bug to the authors"
               "\n\n");
    }

    double log_sum_w = (log_sum_abs_dw_p
                        + log1p(-exp(log_sum_abs_dw_n - log_sum_abs_dw_p))
                        + logscale);
    assert(log_sum_abs_dw_p >= log_sum_abs_dw_n);

    iproc_vector_shift(&logprobs_nz.vector, -log_sum_w);
    *plogsumweight_diff = log_sum_w;
}

static void
iproc_model_ctx_free (iproc_model_ctx *ctx)
{
    if (ctx) {
        iproc_model_unref(ctx->model);
        iproc_svector_unref(ctx->active_probs);
        iproc_svector_unref(ctx->active_logprobs);
        iproc_free(ctx);
    }
}

iproc_model_ctx *
iproc_model_ctx_new (iproc_model     *model,
                     int64_t          isend,
                     iproc_history   *h)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_model_ctx *ctx = iproc_malloc(sizeof(*ctx));

    if (!ctx)
        return NULL;

    ctx->active_probs = NULL;
    ctx->active_logprobs = NULL;

    iproc_vars *vars = iproc_model_vars(model);
    iproc_vars_ctx *vars_ctx = iproc_vars_ctx_new(vars, isend, h);
    iproc_vector *coefs = iproc_model_coefs(model);
    int has_loops = iproc_model_has_loops(model);
    iproc_group_model *group = iproc_model_send_group(model, isend);
    int64_t nreceiver = iproc_model_nreceiver(model);

    ctx->model = iproc_model_ref(model);
    ctx->group = group;
    ctx->active_logprobs = iproc_svector_new(nreceiver);
    
    if (!(vars_ctx && ctx->active_logprobs)) {
        iproc_model_ctx_free(ctx);
        return NULL;
    }

    compute_new_logprobs(vars_ctx, coefs, has_loops, group->logprobs0,
                         ctx->active_logprobs, &ctx->logsumweight_diff);
    iproc_vars_ctx_unref(vars_ctx);

    ctx->active_probs = iproc_svector_new_copy(ctx->active_logprobs);

    if (!ctx->active_probs) {
        iproc_model_ctx_free(ctx);
        return NULL;
    }

    iproc_vector_view view = iproc_svector_view_nz(ctx->active_probs);
    iproc_vector_exp(&view.vector);
    ctx->invsumweight_ratio = exp(-ctx->logsumweight_diff);

    iproc_refcount_init(&ctx->refcount);
    return ctx;
}

iproc_model_ctx *
iproc_model_ctx_ref (iproc_model_ctx *ctx)
{
    if (ctx) {
        iproc_refcount_get(&ctx->refcount);
    }
    return ctx;
}

static void
iproc_model_ctx_release (iproc_refcount *refcount)
{
    iproc_model_ctx *ctx = container_of(refcount, iproc_model_ctx, refcount);
    iproc_model_ctx_free(ctx);
}

void
iproc_model_ctx_unref (iproc_model_ctx *ctx)
{
    if (ctx) {
        iproc_refcount_put(&ctx->refcount, iproc_model_ctx_release);
    }
}

int64_t
iproc_model_ctx_nreceiver (iproc_model_ctx *ctx)
{
    assert(ctx);
    iproc_model *model = ctx->model;
    int64_t nreceiver = iproc_model_nreceiver(model);
    return nreceiver;
}

double
iproc_model_ctx_prob (iproc_model_ctx *ctx,
                      int64_t          jrecv)
{
    assert(ctx);
    assert(0 <= jrecv);
    assert(jrecv < iproc_model_ctx_nreceiver(ctx));

    double lp = iproc_model_ctx_logprob(ctx, jrecv);
    double p = exp(lp);
    return p;
}

double
iproc_model_ctx_logprob (iproc_model_ctx *ctx,
                         int64_t          jrecv)
{
    assert(ctx);
    assert(0 <= jrecv);
    assert(jrecv < iproc_model_ctx_nreceiver(ctx));

    iproc_svector *logprobs = ctx->active_logprobs;
    int64_t jnz = iproc_svector_find_nz(logprobs, jrecv);
    double lp;
    
    if (jnz >= 0) {
        lp = iproc_svector_nz_val(logprobs, jnz);
    } else {
        iproc_vector *logprobs0 = ctx->group->logprobs0;
        double logsumweight_diff = ctx->logsumweight_diff;
        
        lp = iproc_vector_get(logprobs0, jrecv) - logsumweight_diff;
    }

    lp = IPROC_MIN(lp, 0.0);

    return lp;
}

void
iproc_model_ctx_get_probs (iproc_model_ctx *ctx,
                           iproc_vector    *probs)
{
    assert(ctx);
    assert(probs);
    assert(iproc_model_ctx_nreceiver(ctx) == iproc_vector_dim(probs));

    iproc_model_ctx_get_logprobs(ctx, probs);
    iproc_vector_exp(probs);
}

void
iproc_model_ctx_get_logprobs (iproc_model_ctx *ctx,
                              iproc_vector    *logprobs)
{
    assert(ctx);
    assert(logprobs);
    assert(iproc_model_ctx_nreceiver(ctx) == iproc_vector_dim(logprobs));

    int64_t j, n = iproc_model_ctx_nreceiver(ctx);

    for (j = 0; j < n; j++) {
        double lp = iproc_model_ctx_logprob(ctx, j);
        iproc_vector_set(logprobs, j, lp);
    }
}

double
iproc_model_invsumweight (iproc_model_ctx *ctx)
{
    assert(ctx);
    return exp(iproc_model_logsumweight(ctx));
}

double
iproc_model_logsumweight (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->group->logsumweight0 + ctx->logsumweight_diff;
}

double
iproc_model_invsumweight_ratio (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->invsumweight_ratio;
}

double
iproc_model_logsumweight_diff (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->logsumweight_diff;
}

iproc_svector *
iproc_model_active_probs (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->active_probs;
}

iproc_svector *
iproc_model_active_logprobs (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->active_logprobs;
}
