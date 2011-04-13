#include "port.h"

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "model.h"
#include "util.h"

static void
compute_logprobs0(iproc_design * design,
		  int64_t isend,
		  struct vector *coefs,
		  struct vector *logprobs, double *logsumweight)
{
	assert(design);
	assert(coefs);
	assert(logprobs);
	assert(logsumweight);

	/* NOTE: should we worry about overflow in the multiplication?  It is possible
	 * to gaurd against said overflow by scaling the coefficient vector before the
	 * multiplication, then unscaling after subtracting of the max value.  This
	 * shouldn't be necessary in most (all?) real-world situations.
	 */
	iproc_design_mul0(1.0, IPROC_TRANS_NOTRANS, design, isend, coefs,
			  0.0, logprobs);

	/* protect against overflow */
	double max = vector_max(logprobs);
	vector_shift(logprobs, -max);

	/* The probabilities are p[i] = w[i] / sum(w[j]), so
	 * log(p[i]) = log(w[i]) - log(sum(w[j])).
	 */
	double log_sum_exp = vector_log_sum_exp(logprobs);
	vector_shift(logprobs, -log_sum_exp);
	*logsumweight = log_sum_exp + max;
}

static void
iproc_group_models_init(struct darray *group_models,
			iproc_design * design, struct vector *coefs)
{
	int64_t i, nsender = iproc_design_nsender(design);
	int64_t nreceiver = iproc_design_nreceiver(design);
	int64_t dim = iproc_design_dim(design);
	iproc_actors *senders = iproc_design_senders(design);

	darray_init(group_models, iproc_group_model);
	darray_resize(group_models, iproc_actors_ngroup(senders));

	/* We have to loop over all senders, not over all groups, because we
	 * need a representative sender for each group when we call
	 * iproc_design_ctx_new(design, NULL, i).
	 */
	for (i = 0; i < nsender; i++) {
		/* find the group of the sender */
		int64_t g = iproc_actors_group(senders, i);
		iproc_group_model *group = &(darray_index(group_models,
							  iproc_group_model,
							  g));

		/* continue if the group is already initialized */
		if (group->p0)
			continue;

		/* compute initial log(probs) */
		group->log_p0 = vector_new(nreceiver);
		compute_logprobs0(design, i, coefs, group->log_p0,
				  &group->log_W0);

		/* compute initial probs */
		group->p0 = vector_new_copy(group->log_p0);
		vector_exp(group->p0);

		/* compute initial covariate mean */
		group->xbar0 = vector_new(dim);
		iproc_design_mul0(1.0, IPROC_TRANS_TRANS, design, i, group->p0,
				  0.0, group->xbar0);
	}
}

static void iproc_group_models_deinit(struct darray *group_models)
{
	if (!group_models)
		return;

	int64_t i, n = darray_size(group_models);

	for (i = 0; i < n; i++) {
		iproc_group_model *group = &(darray_index(group_models,
							  iproc_group_model,
							  i));
		vector_free(group->xbar0);
		vector_free(group->log_p0);
		vector_free(group->p0);
	}

	darray_deinit(group_models);
}

iproc_group_model *iproc_model_send_group(iproc_model * model, int64_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	struct darray *group_models = &model->group_models;
	iproc_design *design = iproc_model_design(model);
	iproc_actors *senders = iproc_design_senders(design);
	int64_t g = iproc_actors_group(senders, isend);
	return &(darray_index(group_models, iproc_group_model, g));
}

static void iproc_model_ctx_free_dealloc(iproc_model_ctx * ctx)
{
	if (ctx) {
		iproc_svector_unref(ctx->dxbar);
		iproc_svector_unref(ctx->dp);
		iproc_svector_unref(ctx->deta);
		iproc_free(ctx);
	}
}

static void iproc_model_free(iproc_model * model)
{
	if (model) {
		int64_t i, n = darray_size(&model->ctxs);

		for (i = 0; i < n; i++) {
			iproc_model_ctx *ctx = darray_index(&model->ctxs,
							    iproc_model_ctx *,
							    i);
			iproc_model_ctx_free_dealloc(ctx);
		}
		darray_deinit(&model->ctxs);
		iproc_group_models_deinit(&model->group_models);
		vector_free(model->coefs);
		iproc_design_unref(model->design);
		iproc_free(model);
	}
}

iproc_model *iproc_model_new(iproc_design * design,
			     struct vector *coefs, bool has_loops)
{
	assert(design);
	assert(coefs);
	assert(iproc_design_dim(design) == vector_size(coefs));
	assert(iproc_design_nreceiver(design) > 0);
	assert(!has_loops || iproc_design_nreceiver(design) > 1);

	iproc_model *model = iproc_malloc(sizeof(*model));
	model->design = iproc_design_ref(design);
	model->coefs = vector_new_copy(coefs);
	model->has_loops = has_loops;
	iproc_group_models_init(&model->group_models, design, coefs);
	darray_init(&model->ctxs, iproc_model_ctx *);
	refcount_init(&model->refcount);

	return model;
}

iproc_model *iproc_model_ref(iproc_model * model)
{
	if (model) {
		refcount_get(&model->refcount);
	}
	return model;
}

static void iproc_model_release(struct refcount *refcount)
{
	iproc_model *model = container_of(refcount, iproc_model, refcount);
	iproc_model_free(model);
}

void iproc_model_unref(iproc_model * model)
{
	if (!model)
		return;

	refcount_put(&model->refcount, iproc_model_release);
}

iproc_design *iproc_model_design(iproc_model * model)
{
	assert(model);
	return model->design;
}

struct vector *iproc_model_coefs(iproc_model * model)
{
	assert(model);
	return model->coefs;
}

bool iproc_model_has_loops(iproc_model * model)
{
	assert(model);
	return model->has_loops;
}

int64_t iproc_model_nsender(iproc_model * model)
{
	assert(model);
	iproc_design *design = iproc_model_design(model);
	return iproc_design_nsender(design);
}

int64_t iproc_model_nreceiver(iproc_model * model)
{
	assert(model);
	iproc_design *design = iproc_model_design(model);
	return iproc_design_nreceiver(design);
}

int64_t iproc_model_dim(iproc_model * model)
{
	assert(model);
	iproc_design *design = iproc_model_design(model);
	return iproc_design_dim(design);
}

struct vector *iproc_model_logprobs0(iproc_model * model, int64_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_group_model *group = iproc_model_send_group(model, isend);
	assert(group);
	return group->log_p0;
}

struct vector *iproc_model_probs0(iproc_model * model, int64_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_group_model *group = iproc_model_send_group(model, isend);
	assert(group);
	return group->p0;
}

struct vector *iproc_model_mean0(iproc_model * model, int64_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_group_model *group = iproc_model_send_group(model, isend);
	assert(group);
	return group->xbar0;
}
