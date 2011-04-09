#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "logsumexp.h"
#include "memory.h"
#include "refcount.h"
#include "model.h"


static void
compute_weight_changes (iproc_design_ctx *ctx,
                        bool              has_loops,
                        struct vector     *coefs,
                        struct vector     *log_p0,
                        iproc_svector    *deta,
                        double           *gamma,
                        double           *log_gamma)
{
    /* compute the changes in weights */
    iproc_design_ctx_dmul(1.0, IPROC_TRANS_NOTRANS, ctx, coefs, 0.0, deta);
    if (!has_loops) {
        iproc_svector_set(deta, ctx->isend, -INFINITY);
    }
    
    int64_t i, nnz = iproc_svector_nnz(deta);
    iproc_vector_view deta_nz = iproc_svector_view_nz(deta);
    
    /* compute the scale for the weight differences */
    double lwmax = vector_max(&deta_nz.vector);
    double logscale = MAX(0.0, lwmax);
    double invscale = exp(-logscale);
    
    iproc_logsumexp pos, neg;
    iproc_logsumexp_init(&pos);
    iproc_logsumexp_init(&neg);
    
    /* treat the initial sum of weights as the first positive difference */
    iproc_logsumexp_insert(&pos, -logscale);
    
    /* compute the log sums of the positive and negative differences in weights */
    for (i = 0; i < nnz; i++) {
        int64_t jrecv = iproc_svector_nz(deta, i);
        double lp0 = vector_index(log_p0, jrecv);
        double dlw = vector_index(&deta_nz.vector, i);
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
    
    *gamma = exp(-log_sum_w);
    *log_gamma = -log_sum_w;
}

/* Given log_p0, log_gamma, and deta, compute p_active.  The probability for
 * receiver j is active if x[t,i,j] is nonzero or j = i and self-loops are
 * dis-allowed.  We rely on the fact that deta and p_active have the same
 * sparsity pattern.
 */
static void
compute_active_probs (struct vector  *log_p0,
                      double         log_gamma,
                      iproc_svector *deta,
                      iproc_svector *p_active)
{
    iproc_svector_copy(p_active, deta);
    int64_t i, nnz = iproc_svector_nnz(p_active);
    iproc_vector_view p_active_nz = iproc_svector_view_nz(p_active);
    
    for (i = 0; i < nnz; i++) {
        int64_t jrecv = iproc_svector_nz(p_active, i);
        double lp0_j = vector_index(log_p0, jrecv);
        double deta_j = vector_index(&p_active_nz.vector, i);
        double lp_j = MIN(0.0, log_gamma + lp0_j + deta_j);
        vector_index(&p_active_nz.vector, i) = lp_j;
    }
    
    vector_exp(&p_active_nz.vector);
}

static void
compute_prob_diffs(struct vector  *p0,
                   double         gamma,
                   iproc_svector *p_active)
{
    int64_t i, nnz = iproc_svector_nnz(p_active);
    iproc_vector_view p_active_nz = iproc_svector_view_nz(p_active);
    
    for (i = 0; i < nnz; i++) {
        int64_t jrecv = iproc_svector_nz(p_active, i);
        double p0_j = vector_index(p0, jrecv);
        double p_j = vector_index(&p_active_nz.vector, i);
        double dp_j = p_j - gamma * p0_j;
        vector_index(&p_active_nz.vector, i) = dp_j;
    }
}

static void
iproc_model_ctx_free (iproc_model_ctx *ctx)
{
    if (ctx) {
        iproc_design_ctx_unref(ctx->design_ctx);

        iproc_model *model = ctx->model;
        darray_push_back(&model->ctxs, &ctx);
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
    int64_t dim = iproc_model_dim(model);

    ctx->model = iproc_model_ref(model);
    ctx->design_ctx = NULL;
    ctx->group = NULL;
    
    ctx->deta = iproc_svector_new(nreceiver);
    ctx->dp = iproc_svector_new(nreceiver);
    ctx->dxbar = iproc_svector_new(dim);
    
    iproc_refcount_init(&ctx->refcount);
    
    if (!(ctx->deta && ctx->dp && ctx->dxbar)) {
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
    struct darray *ctxs = &model->ctxs;
    int64_t n = darray_size(ctxs);
    
    if (n > 0) {
        ctx = darray_index(ctxs, iproc_model_ctx *, n - 1);
        darray_resize(ctxs, n - 1);
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

    iproc_svector_clear(ctx->deta);
    iproc_svector_clear(ctx->dp);
    iproc_svector_clear(ctx->dxbar);
    iproc_design_ctx_unref(ctx->design_ctx);

    iproc_model *model = ctx->model;
    iproc_design *design = iproc_model_design(model);
    iproc_design_ctx *design_ctx = iproc_design_ctx_new(design, isend, h);
    struct vector *coefs = iproc_model_coefs(model);
    bool has_loops = iproc_model_has_loops(model);
    iproc_group_model *group = iproc_model_send_group(model, isend);

    ctx->design_ctx = design_ctx;
    ctx->group = group;

    compute_weight_changes(design_ctx, has_loops, coefs, group->log_p0,
                           ctx->deta, &ctx->gamma, &ctx->log_gamma);
    compute_active_probs(group->log_p0, ctx->log_gamma, ctx->deta, ctx->dp);
    iproc_design_ctx_dmuls(1.0, IPROC_TRANS_TRANS, design_ctx, ctx->dp,
                               0.0, ctx->dxbar);
    compute_prob_diffs(group->p0, ctx->gamma, ctx->dp);
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

/*
 *         log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double
iproc_model_ctx_logprob (iproc_model_ctx *ctx,
                         int64_t          jrecv)
{
    assert(ctx);
    assert(0 <= jrecv);
    assert(jrecv < iproc_model_ctx_nreceiver(ctx));

    /*
    double gamma = ctx->gamma;
    double p0 = iproc_vector_get(ctx->group->p0, jrecv);
    double dp = iproc_svector_get(ctx->dp, jrecv);
    double p = gamma * p0 + dp;
    p = IPROC_MAX(0.0, IPROC_MIN(1.0, p));
    return log(p);
     */
    
    double log_gamma = ctx->log_gamma;
    double log_p0 = vector_index(ctx->group->log_p0, jrecv);
    double deta = iproc_svector_get(ctx->deta, jrecv);
    double log_p = log_gamma + log_p0 + deta;
    return MIN(log_p, 0.0);
}

void
iproc_model_ctx_get_probs (iproc_model_ctx *ctx,
                           struct vector    *probs)
{
    assert(ctx);
    assert(probs);
    assert(iproc_model_ctx_nreceiver(ctx) == vector_size(probs));

    iproc_model_ctx_get_logprobs(ctx, probs);
    vector_exp(probs);
}

void
iproc_model_ctx_get_logprobs (iproc_model_ctx *ctx,
                              struct vector    *logprobs)
{
    assert(ctx);
    assert(logprobs);
    assert(iproc_model_ctx_nreceiver(ctx) == vector_size(logprobs));

    int64_t j, n = iproc_model_ctx_nreceiver(ctx);

    for (j = 0; j < n; j++) {
        double lp = iproc_model_ctx_logprob(ctx, j);
        vector_index(logprobs, j) = lp;
    }
}
