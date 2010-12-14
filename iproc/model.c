#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include <iproc/memory.h>
#include <iproc/model.h>


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
