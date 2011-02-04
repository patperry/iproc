#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "model.h"
#include "utils.h"

static void
compute_logprobs0 (iproc_frame   *frame,
                   int64_t       isend,
                   iproc_vector *coefs,
                   int           has_loops,
                   iproc_vector *logprobs,
                   double       *logsumweight)
{
    assert(frame);
    assert(coefs);
    assert(logprobs);
    assert(logsumweight);

    /* NOTE: should we worry about overflow in the multiplication?  It is possible
     * to gaurd against said overflow by scaling the coefficient vector before the
     * multiplication, then unscaling after subtracting of the max value.  This
     * shouldn't be necessary in most (all?) real-world situations.
     */
    iproc_frame_sender0_mul(1.0, IPROC_TRANS_NOTRANS, frame, isend, coefs,
                           0.0, logprobs);

    if (!has_loops) {
        iproc_vector_set(logprobs, isend, -INFINITY);
    }

    /* protect against overflow */
    double max = iproc_vector_max(logprobs);
    iproc_vector_shift(logprobs, -max);

    /* The probabilities are p[i] = w[i] / sum(w[j]), so
     * log(p[i]) = log(w[i]) - log(sum(w[j])).
     */
    double log_sum_exp = iproc_vector_log_sum_exp(logprobs);
    iproc_vector_shift(logprobs, -log_sum_exp);
    *logsumweight = log_sum_exp + max;
}


static void
iproc_group_models_init (iproc_array  *group_models,
                         iproc_frame   *frame,
                         iproc_vector *coefs)
{
    int64_t i, nsender = iproc_frame_nsender(frame);
    int64_t nreceiver = iproc_frame_nreceiver(frame);
    int64_t dim = iproc_frame_dim(frame);
    iproc_actors *senders = iproc_frame_senders(frame);

    iproc_array_set_size(group_models, iproc_actors_ngroup(senders));

    /* We have to loop over all senders, not over all groups, because we
     * need a representative sender for each group when we call
     * iproc_frame_ctx_new(frame, NULL, i).
     */
    for (i = 0; i < nsender; i++) {
        int64_t k = iproc_actors_group(senders, i);
        iproc_group_model *group = &(iproc_array_index(group_models,
                                                       iproc_group_model,
                                                       k));
        if (group->logprobs0)
            continue;

        group->logprobs0 = iproc_vector_new(nreceiver);
        compute_logprobs0(frame, i, coefs, 1, group->logprobs0, &group->logsumweight0);

        group->probs0 = iproc_vector_new_copy(group->logprobs0);
        iproc_vector_exp(group->probs0);
        group->invsumweight0 = exp(-group->logsumweight0);

        group->mean0 = iproc_vector_new(dim);
        iproc_frame_sender0_mul(1.0, IPROC_TRANS_TRANS, frame, i, group->probs0,
                               0.0, group->mean0);
    }
}


static void
iproc_group_models_deinit (iproc_array *group_models)
{
    if (!group_models)
        return;

    int64_t n = iproc_array_size(group_models);
    int64_t i;

    for (i = 0; i < n; i++) {
        iproc_group_model *group = &(iproc_array_index(group_models,
                                                       iproc_group_model,
                                                       i));
        iproc_vector_unref(group->logprobs0);
        iproc_vector_unref(group->probs0);
        iproc_vector_unref(group->mean0);
    }
}

iproc_group_model *
iproc_model_send_group (iproc_model *model,
                        int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_array *group_models = model->group_models;
    iproc_frame *frame = iproc_model_frame(model);
    iproc_actors *senders = iproc_frame_senders(frame);
    int64_t g = iproc_actors_group(senders, isend);
    return &(iproc_array_index(group_models, iproc_group_model, g));
}

static void
iproc_model_free (iproc_model *model)
{
    if (model) {
        iproc_vector_unref(model->coefs);
        iproc_frame_unref(model->frame);
        iproc_group_models_deinit(model->group_models);
        iproc_array_unref(model->group_models);
        iproc_free(model);
    }
}

iproc_model *
iproc_model_new (iproc_frame   *frame,
                 iproc_vector *coefs,
                 int           has_loops)
{
    assert(frame);
    assert(coefs);
    assert(iproc_frame_dim(frame) == iproc_vector_dim(coefs));
    assert(iproc_frame_nreceiver(frame) > 0);
    assert(!has_loops || iproc_frame_nreceiver(frame) > 1);

    iproc_model *model = iproc_malloc(sizeof(*model));
    model->frame = iproc_frame_ref(frame);
    model->coefs = iproc_vector_new_copy(coefs);
    model->has_loops = has_loops;
    model->group_models = iproc_array_new(sizeof(iproc_group_model));
    iproc_group_models_init(model->group_models, frame, coefs);
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

iproc_frame *
iproc_model_frame (iproc_model *model)
{
    assert(model);
    return model->frame;
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
    iproc_frame *frame = iproc_model_frame(model);
    return iproc_frame_nsender(frame);
}

int64_t
iproc_model_nreceiver (iproc_model *model)
{
    assert(model);
    iproc_frame *frame = iproc_model_frame(model);
    return iproc_frame_nreceiver(frame);
}

int64_t
iproc_model_dim (iproc_model *model)
{
    assert(model);
    iproc_frame *frame = iproc_model_frame(model);
    return iproc_frame_dim(frame);
}

double
iproc_model_invsumweight0 (iproc_model *model,
                           int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_group_model *group = iproc_model_send_group(model, isend);
    assert(group);

    return group->invsumweight0;
}

double
iproc_model_logsumweight0 (iproc_model *model,
                           int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_group_model *group = iproc_model_send_group(model, isend);
    assert(group);

    return group->logsumweight0;
}

iproc_vector *
iproc_model_logprobs0 (iproc_model *model,
                       int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_group_model *group = iproc_model_send_group(model, isend);
    assert(group);
    return group->logprobs0;
}

iproc_vector *
iproc_model_probs0 (iproc_model *model,
                    int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_group_model *group = iproc_model_send_group(model, isend);
    assert(group);
    return group->probs0;
}

iproc_vector *
iproc_model_mean0 (iproc_model *model,
                    int64_t      isend)
{
    assert(model);
    assert(isend >= 0);
    assert(isend < iproc_model_nsender(model));

    iproc_group_model *group = iproc_model_send_group(model, isend);
    assert(group);
    return group->mean0;
}
