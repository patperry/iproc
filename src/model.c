#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "model.h"
#include "utils.h"

static iproc_array *
iproc_model_array_new (iproc_vars   *vars,
                       iproc_vector *coef)
{
    assert(vars);
    assert(coef);

    iproc_actors *senders = iproc_vars_senders(vars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    int64_t i, nsender = iproc_actors_size(senders);
    int64_t nreceiver = iproc_actors_size(receivers);

    iproc_array *array = iproc_array_new(sizeof(iproc_vector *));
    iproc_array_set_size(array, iproc_actors_ngroup(senders));

    for (i = 0; i < nsender; i++) {
        int64_t k = iproc_actors_group(senders, i);
        iproc_vector *logprobs = iproc_array_index(array, iproc_vector *, k);

        if (logprobs)
            continue;

        iproc_vars_ctx *ctx = iproc_vars_ctx_new(vars, NULL, i);

        logprobs = iproc_vector_new(nreceiver);
        iproc_vars_ctx_mul(1.0, IPROC_TRANS_NOTRANS, ctx, coef, 0.0, logprobs);
        iproc_vector_shift(logprobs, -iproc_vector_max(logprobs));
        double log_sum_exp = iproc_vector_log_sum_exp(logprobs);
        iproc_vector_shift(logprobs, -log_sum_exp);

        iproc_array_index(array, iproc_vector *, k) = logprobs;
        iproc_vars_ctx_unref(ctx);
    }

    return array;
}

static void
iproc_model_array_unref (iproc_array *array)
{
    if (!array)
        return;

    int64_t n = iproc_array_size(array);
    int64_t i;

    for (i = 0; i < n; i++) {
        iproc_vector *vector = iproc_array_index(array, iproc_vector *, i);
        iproc_vector_unref(vector);
    }

    iproc_array_unref(array);
}

static void
iproc_model_free (iproc_model *model)
{
    if (model) {
        iproc_vector_unref(model->coef);
        iproc_vars_unref(model->vars);
        iproc_model_array_unref(model->logprobs0_array);
        iproc_free(model);
    }
}

iproc_model *
iproc_model_new (iproc_vars   *vars,
                 iproc_vector *coef,
                 int           has_loops)
{
    assert(vars);
    assert(coef);
    assert(iproc_vars_dim(vars) == iproc_vector_dim(coef));
    assert(iproc_vars_nreceiver(vars) > 0);
    assert(!has_loops || iproc_vars_nreceiver(vars) > 1);

    iproc_model *model = iproc_malloc(sizeof(*model));
    model->vars = iproc_vars_ref(vars);
    model->coef = iproc_vector_new_copy(coef);
    model->has_loops = has_loops;
    iproc_refcount_init(&model->refcount);
    model->logprobs0_array = iproc_model_array_new(vars, coef);
    return model;
}

iproc_model *
iproc_model_ref (iproc_model *model)
{
    if (model) {
        iproc_refcount_get(&model->refcount);
    }
    return model;
}

static void
iproc_model_release (iproc_refcount *refcount)
{
    iproc_model *model = container_of(refcount, iproc_model, refcount);
    iproc_model_free(model);
}

void
iproc_model_unref (iproc_model *model)
{
    if (!model)
        return;

    iproc_refcount_put(&model->refcount, iproc_model_release);
}

void
iproc_model_get_probs (iproc_model    *model,
                       iproc_vars_ctx *ctx,
                       iproc_vector   *probs)
{
    iproc_model_get_logprobs(model, ctx, probs);
    iproc_vector_exp(probs);
}

void
iproc_model_get_logprobs (iproc_model    *model,
                          iproc_vars_ctx *ctx,
                          iproc_vector   *logprobs)
{
    assert(model);
    assert(ctx);
    assert(logprobs);
    assert(model->vars == ctx->vars);
    assert(iproc_vars_nreceiver(model->vars) == iproc_vector_dim(logprobs));

    int64_t imax, i, n = iproc_vector_dim(logprobs);
    double max, summ1;

    /* NOTE: should we worry about overflow in the multiplication?  It is possible
     * to gaurd against said overflow by scaling the coefficient vector before the
     * multiplication, then unscaling after subtracting of the max value.  This
     * shouldn't be necessary in most (all?) real-world situations.
     */
    iproc_vars_ctx_mul(1.0, IPROC_TRANS_NOTRANS, ctx, model->coef, 0.0, logprobs);

    if (!model->has_loops) {
        iproc_vector_set(logprobs, ctx->isend, -INFINITY);
    }

    imax = iproc_vector_max_index(logprobs);
    max = iproc_vector_get(logprobs, imax);
    iproc_vector_shift(logprobs, -max);

    /* compute the (scaled) sum of the weights, minus 1 */
    summ1 = 0.0;
    for (i = 0; i < n; i++) {
        if (i == imax)
            continue;
        summ1 += exp(iproc_vector_get(logprobs, i));
    }

    /* The probabilities are p[i] = w[i] / sum(w[j]), so
     * log(p[i]) = log(w[i]) - log(sum(w[j])).
     */
    iproc_vector_shift(logprobs, -log1p(summ1));
}


double
iproc_model_logprob0 (iproc_model *model,
                      int64_t      i,
                      int64_t      j)
{
    iproc_vars *vars = model->vars;
    iproc_actors *senders = iproc_vars_senders(vars);
    int64_t k = iproc_actors_group(senders, i);
    iproc_array *array = model->logprobs0_array;
    iproc_vector *logprobs0 = iproc_array_index(array, iproc_vector *, k);
    double lp0 = iproc_vector_get(logprobs0, j);
    return lp0;
}

/* Algorithm logsumexp
 *   INPUT: a[1], ..., a[N]
 *   OUTPUT: log(sum(exp(ai)))
 *
 *   INVARIANT: S[n] = sum(exp(a[i] - M[n])) - 1
 *              M[n] = max(a[i])
 *
 *   S := -1
 *   M := -INFINITY
 *
 *   for i = 1 to N do
 *       if a[i] <= M
 *       then
 *           S := S + exp(a[i] - M)
 *       else
 *           S := (S + 1) * exp(M - a[i])
 *           M := a[i]
 *       end
 *   end
 *
 *   return M + log(S + 1)
 */
static void
logsumexp_init (double *psum_exp, double *pmax)
{
    *psum_exp = -1.0;
    *pmax = -INFINITY;
}


static void
logsumexp_iter (double val, double *psum_exp, double *pmax)
{
    double sum_exp = *psum_exp;
    double max = *pmax;

    if (val <= *pmax) {
        *psum_exp = sum_exp + exp(val - max);
    } else {
        double scale = exp(max - val);
        *psum_exp = sum_exp * scale + scale;
        *pmax = val;
    }
}

static double
logsumexp_result (double *psum_exp, double *pmax)
{
    return (*pmax + log1p(*psum_exp));
}

void
iproc_model_get_new_logprobs (iproc_model    *model,
                              iproc_vars_ctx *ctx,
                              double         *plogprob0_shift,
                              iproc_svector  *logprobs)
{
    assert(model);
    assert(ctx);
    assert(logprobs);
    assert(model->vars == ctx->vars);
    assert(plogprob0_shift);
    assert(iproc_vars_nreceiver(model->vars) == iproc_svector_dim(logprobs));

    int64_t nnz, i;

    /* compute the changes in weights */
    iproc_vars_ctx_diff_mul(1.0, IPROC_TRANS_NOTRANS, ctx, model->coef, 0.0, logprobs);
    nnz = iproc_svector_nnz(logprobs);

    if (nnz == 0) {
        *plogprob0_shift = 0;
        return;
    }

    iproc_vector_view logprobs_nz = iproc_svector_view_nz(logprobs);

    /* compute the scale for the weight differences */
    double lwmax = iproc_vector_max(&logprobs_nz.vector);
    double logscale = IPROC_MAX(0.0, lwmax);
    double invscale = exp(-logscale);

    double sp, mp, sn, mn;
    logsumexp_init(&sp, &mp);
    logsumexp_init(&sn, &mn);

    /* treat the initial sum of weights as the first positive difference */
    logsumexp_iter(-logscale, &sp, &mp);

    /* treat the self-loop as the first negative difference */
    if (!model->has_loops) {
        double lp0 = iproc_model_logprob0(model, i, i);
        double log_abs_dw = lp0 + invscale;
        logsumexp_iter(log_abs_dw, &sn, &mn);
    }

    /* compute the log sums of the positive and negative differences in weights */
    for (i = 0; i < nnz; i++) {
        int64_t j = iproc_svector_nz(logprobs, i);
        double lp0 = iproc_model_logprob0(model, i, j);
        double dlw = iproc_vector_get(&logprobs_nz.vector, i);
        double log_abs_dw;

        /* special handling for self-loops */
        if (j == i && !model->has_loops) {
            /* no need to call logsumexp_iter since self-loop gets handled before
             * the for loop */
            iproc_vector_set(&logprobs_nz.vector, i, -INFINITY);
            continue;
        }

        if (dlw >= 0) {
            log_abs_dw = lp0 + (exp(dlw - logscale) - invscale);
            logsumexp_iter(log_abs_dw, &sp, &mp);
        } else {
            log_abs_dw = lp0 + (invscale - exp(dlw - logscale));
            logsumexp_iter(log_abs_dw, &sn, &mn);
        }

        /* make logprobs[j] store the (unnormalized) log probability */
        iproc_vector_inc(&logprobs_nz.vector, i, lp0);
    }

    double log_sum_abs_dw_p = logsumexp_result(&sp, &mp);
    double log_sum_abs_dw_n = logsumexp_result(&sn, &mn);
    double log_sum_w = (log_sum_abs_dw_p
                        + log1p(-exp(log_sum_abs_dw_n - log_sum_abs_dw_p))
                        + logscale);
    assert(log_sum_abs_dw_p >= log_sum_abs_dw_n);

    iproc_vector_shift(&logprobs_nz.vector, -log_sum_w);
    *plogprob0_shift = -log_sum_w;
}