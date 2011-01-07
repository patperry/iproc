#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include <iproc/memory.h>
#include <iproc/model.h>
#include <iproc/utils.h>


static void
iproc_model_free (iproc_model *model)
{
    if (model) {
        iproc_vector_unref(model->coef);
        iproc_vars_unref(model->vars);
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
    assert(iproc_vars_nrecv(vars) > 0);
    assert(!has_loops || iproc_vars_nrecv(vars) > 1);

    iproc_model *model = iproc_malloc(sizeof(*model));
    model->vars = iproc_vars_ref(vars);
    model->coef = iproc_vector_new_copy(coef);
    model->has_loops = has_loops;
    model->refcount = 1;
    return model;
}

iproc_model *
iproc_model_ref (iproc_model *model)
{
    if (model) {
        model->refcount = model->refcount + 1;
    }
    return model;
}

void
iproc_model_unref (iproc_model *model)
{
    if (!model)
        return;

    if (model->refcount == 1) {
        iproc_model_free(model);
    } else {
        model->refcount = model->refcount - 1;
    }
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
    assert(iproc_vars_nrecv(model->vars) == iproc_vector_dim(logprobs));

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
    assert(iproc_vars_nrecv(model->vars) == iproc_svector_dim(logprobs));

    int64_t nnz, i;

    /* compute the changes in weights */
    iproc_vars_ctx_diff_mul(1.0, IPROC_TRANS_NOTRANS, ctx, model->coef, 0.0, logprobs);
    nnz = iproc_svector_nnz(logprobs);

    if (nnz == 0) {
        *plogprob0_shift = 0;
        return;
    }

    /* compute the scale for the weight differences */
    double lwmax = iproc_svector_nz_max(logprobs);
    double logscale = IPROC_MAX(0.0, lwmax);
    double invscale = exp(-logscale);

    double sp, mp, sn, mn;
    logsumexp_init(&sp, &mp);
    logsumexp_init(&sn, &mn);

    /* treat the initial sum of weights as the first positive difference */
    logsumexp_iter(-logscale, &sp, &mp);

    /* compute the log sums of the positive and negative differences in weights */
    for (i = 0; i < nnz; i++) {
        int64_t j = iproc_svector_nz(logprobs, i);
        double lp0 = iproc_model_logprob0(model, j);
        double dlw = iproc_svector_nz_val(logprobs, i);
        double log_abs_dw;

        if (dlw >= 0) {
            log_abs_dw = lp0 + (exp(dlw - logscale) - invscale);
            logsumexp_iter(log_abs_dw, &sp, &mp);
        } else {
            log_abs_dw = lp0 + (invscale - exp(dlw - logscale));
            logsumexp_iter(log_abs_dw, &sn, &mn);
        }

        /* make logprobs[j] store the (unnormalized) log probability */
        iproc_svector_nz_inc(logprobs, i, lp0);
    }

    double log_sum_abs_dw_p = logsumexp_result(&sp, &mp);
    double log_sum_abs_dw_n = logsumexp_result(&sn, &mn);
    double log_sum_w = (log_sum_abs_dw_p
                        + log1p(-exp(log_sum_abs_dw_n - log_sum_abs_dw_p))
                        + logscale);
    assert(log_sum_abs_dw_p >= log_sum_abs_dw_n);

    iproc_svector_nz_shift(logprobs, -log_sum_w);
    *plogprob0_shift = -log_sum_w;
}
