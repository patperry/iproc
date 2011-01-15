#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "model.h"
#include "utils.h"

static void
compute_logprobs (iproc_vars_ctx *ctx,
                  iproc_vector   *coefs,
                  int             has_loops,
                  iproc_vector   *logprobs)
{
    assert(ctx);
    assert(coefs);
    assert(logprobs);

    /* NOTE: should we worry about overflow in the multiplication?  It is possible
     * to gaurd against said overflow by scaling the coefficient vector before the
     * multiplication, then unscaling after subtracting of the max value.  This
     * shouldn't be necessary in most (all?) real-world situations.
     */
    iproc_vars_ctx_mul(1.0, IPROC_TRANS_NOTRANS, ctx, coefs, 0.0, logprobs);

    if (!has_loops) {
        iproc_vector_set(logprobs, ctx->isend, -INFINITY);
    }

    /* protect against overflow */
    double max = iproc_vector_max(logprobs);
    iproc_vector_shift(logprobs, -max);

    /* The probabilities are p[i] = w[i] / sum(w[j]), so
     * log(p[i]) = log(w[i]) - log(sum(w[j])).
     */
    double log_sum_exp = iproc_vector_log_sum_exp(logprobs);
    iproc_vector_shift(logprobs, -log_sum_exp);
}


static void
iproc_group_logprobs0_init (iproc_array  *group_logprobs0,
                            iproc_vars   *vars,
                            iproc_vector *coefs)
{
    int64_t i, nsender = iproc_vars_nsender(vars);
    int64_t nreceiver = iproc_vars_nreceiver(vars);
    iproc_actors *senders = iproc_vars_senders(vars);

    iproc_array_set_size(group_logprobs0, iproc_actors_ngroup(senders));

    /* We have to loop over all senders, not over all groups, because we
     * need a representative sender for each group when we call
     * iproc_vars_ctx_new(vars, NULL, i).
     */
    for (i = 0; i < nsender; i++) {
        int64_t k = iproc_actors_group(senders, i);
        iproc_vector *logprobs0 = iproc_array_index(group_logprobs0,
                                                    iproc_vector *,
                                                    k);

        if (logprobs0)
            continue;

        iproc_vars_ctx *ctx = iproc_vars_ctx_new(vars, i, NULL);

        logprobs0 = iproc_vector_new(nreceiver);
        compute_logprobs(ctx, coefs, 1, logprobs0);
        iproc_array_index(group_logprobs0, iproc_vector *, k) = logprobs0;
        iproc_vars_ctx_unref(ctx);
    }
}


static void
iproc_group_logprobs0_deinit (iproc_array *group_logprobs0)
{
    if (!group_logprobs0)
        return;

    int64_t n = iproc_array_size(group_logprobs0);
    int64_t i;

    for (i = 0; i < n; i++) {
        iproc_vector *logprobs0 = iproc_array_index(group_logprobs0, iproc_vector *, i);
        iproc_vector_unref(logprobs0);
    }
}

static void
iproc_model_free (iproc_model *model)
{
    if (model) {
        iproc_vector_unref(model->coefs);
        iproc_vars_unref(model->vars);
        iproc_group_logprobs0_deinit(model->group_logprobs0);
        iproc_array_unref(model->group_logprobs0);
        iproc_free(model);
    }
}

iproc_model *
iproc_model_new (iproc_vars   *vars,
                 iproc_vector *coefs,
                 int           has_loops)
{
    assert(vars);
    assert(coefs);
    assert(iproc_vars_dim(vars) == iproc_vector_dim(coefs));
    assert(iproc_vars_nreceiver(vars) > 0);
    assert(!has_loops || iproc_vars_nreceiver(vars) > 1);

    iproc_model *model = iproc_malloc(sizeof(*model));
    model->vars = iproc_vars_ref(vars);
    model->coefs = iproc_vector_new_copy(coefs);
    model->has_loops = has_loops;
    model->group_logprobs0 = iproc_array_new(sizeof(iproc_vector *));
    iproc_group_logprobs0_init(model->group_logprobs0, vars, coefs);
    iproc_refcount_init(&model->refcount);

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

iproc_vars *
iproc_model_vars (iproc_model *model)
{
    assert(model);
    return model->vars;
}

iproc_vector *
iproc_model_coefs (iproc_model *model)
{
    assert(model);
    return model->coefs;
}

int
iproc_model_has_loops (iproc_model *model)
{
    assert(model);
    return model->has_loops;
}

int64_t
iproc_model_nsender (iproc_model *model)
{
    assert(model);
    iproc_vars *vars = iproc_model_vars(model);
    return iproc_vars_nsender(vars);
}

int64_t
iproc_model_nreceiver (iproc_model *model)
{
    assert(model);
    iproc_vars *vars = iproc_model_vars(model);
    return iproc_vars_nreceiver(vars);
}

int64_t
iproc_model_dim (iproc_model *model)
{
    assert(model);
    iproc_vars *vars = iproc_model_vars(model);
    return iproc_vars_dim(vars);
}

iproc_vector *
iproc_model_logprobs0 (iproc_model *model,
                       int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_array *group_logprobs0 = model->group_logprobs0;
    iproc_vars *vars = iproc_model_vars(model);
    iproc_actors *senders = iproc_vars_senders(vars);
    int64_t g = iproc_actors_group(senders, isend);
    iproc_vector *lp0 = iproc_array_index(group_logprobs0,
                                          iproc_vector *,
                                          g);
    return lp0;
}
