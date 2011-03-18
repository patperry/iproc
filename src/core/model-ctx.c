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
compute_new_logprobs (iproc_design_ctx *ctx,
                      iproc_vector     *coefs,
                      bool              has_loops,
                      iproc_vector     *logprobs0,
                      iproc_svector    *logprobs,
                      double           *plogsumweight_diff)
{
    assert(ctx);
    assert(coefs);
    assert(logprobs0);
    assert(logprobs);
    assert(plogsumweight_diff);
    assert(iproc_vector_dim(logprobs0) == iproc_svector_dim(logprobs));

    /* compute the changes in weights */
    iproc_design_ctx_diff_mul(1.0, IPROC_TRANS_NOTRANS, ctx, coefs, 0.0, logprobs);
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
        iproc_design_ctx_unref(ctx->design_ctx);

        iproc_model *model = ctx->model;
        iproc_array_append(model->ctxs, &ctx);
        iproc_model_unref(model);
    }
}

static iproc_model_ctx *
iproc_model_ctx_new_alloc (iproc_model     *model,
                           int64_t          isend,
                           iproc_history   *h)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_model_ctx *ctx = iproc_malloc(sizeof(*ctx));
    if (!ctx)
        return NULL;

    int64_t nreceiver = iproc_model_nreceiver(model);
    ctx->model = iproc_model_ref(model);
    ctx->design_ctx = NULL;
    ctx->group = NULL;
    ctx->active_logprobs = iproc_svector_new(nreceiver);
    ctx->active_probs = iproc_svector_new(nreceiver);
    iproc_refcount_init(&ctx->refcount);
    
    if (!(ctx->active_logprobs && ctx->active_probs)) {
        iproc_model_ctx_free(ctx);
        return NULL;
    }

    iproc_model_ctx_set(ctx, isend, h);
    return ctx;
}


iproc_model_ctx *
iproc_model_ctx_new (iproc_model     *model,
                     int64_t          isend,
                     iproc_history   *h)
{
    assert(model);
    assert(0 <= isend);
    assert(isend < iproc_model_nsender(model));

    iproc_model_ctx *ctx;
    iproc_array *ctxs = model->ctxs;
    int64_t n = iproc_array_size(ctxs);
    
    if (n > 0) {
        ctx = iproc_array_index(ctxs, iproc_model_ctx *, n - 1);
        iproc_array_set_size(ctxs, n - 1);
        iproc_model_ref(model);
        iproc_refcount_init(&ctx->refcount);
        iproc_model_ctx_set(ctx, isend, h);
        assert(ctx->model == model);
    } else {
        ctx = iproc_model_ctx_new_alloc(model, isend, h);
    }

    return ctx;

}


void
iproc_model_ctx_set (iproc_model_ctx *ctx,
                     int64_t          isend,
                     iproc_history   *h)
{
    assert(ctx);
    assert(ctx->model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(ctx->model));

    iproc_design_ctx_unref(ctx->design_ctx);
    iproc_svector_clear(ctx->active_probs);
    iproc_svector_clear(ctx->active_logprobs);

    iproc_model *model = ctx->model;
    iproc_design *design = iproc_model_design(model);
    iproc_design_ctx *design_ctx = iproc_design_ctx_new(design, isend, h);
    iproc_vector *coefs = iproc_model_coefs(model);
    bool has_loops = iproc_model_has_loops(model);
    iproc_group_model *group = iproc_model_send_group(model, isend);

    ctx->design_ctx = design_ctx;
    ctx->group = group;

    compute_new_logprobs(design_ctx, coefs, has_loops, group->logprobs0,
                         ctx->active_logprobs, &ctx->logsumweight_diff);
    iproc_svector_copy(ctx->active_probs, ctx->active_logprobs);
    
    iproc_vector_view view = iproc_svector_view_nz(ctx->active_probs);
    iproc_vector_exp(&view.vector);
    ctx->invsumweight_ratio = exp(-ctx->logsumweight_diff);
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
        lp = iproc_svector_nz_get(logprobs, jnz);
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
iproc_model_ctx_invsumweight (iproc_model_ctx *ctx)
{
    assert(ctx);
    return exp(iproc_model_ctx_logsumweight(ctx));
}

double
iproc_model_ctx_logsumweight (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->group->logsumweight0 + ctx->logsumweight_diff;
}

double
iproc_model_ctx_invsumweight_ratio (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->invsumweight_ratio;
}

double
iproc_model_ctx_logsumweight_diff (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->logsumweight_diff;
}

iproc_svector *
iproc_model_ctx_active_probs (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->active_probs;
}

iproc_svector *
iproc_model_ctx_active_logprobs (iproc_model_ctx *ctx)
{
    assert(ctx);
    return ctx->active_logprobs;
}
