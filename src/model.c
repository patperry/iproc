#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "model.h"
#include "util.h"

static void
compute_logprobs0(const struct design *design,
		  ssize_t isend,
		  const struct vector *coefs,
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
	design_mul0(1.0, TRANS_NOTRANS, design, isend, coefs, 0.0, logprobs);

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

static bool cohort_model_init(struct cohort_model *cm,
			      const struct design *design,
			      ssize_t isend, const struct vector *coefs)
{
	ssize_t nreceiver = design_nreceiver(design);
	ssize_t dim = design_dim(design);

	/* compute initial log(probs) */
	cm->log_p0 = vector_alloc(nreceiver);
	compute_logprobs0(design, isend, coefs, cm->log_p0, &cm->log_W0);

	/* compute initial probs */
	cm->p0 = vector_alloc_copy(cm->log_p0);
	vector_exp(cm->p0);

	/* compute initial covariate mean */
	cm->xbar0 = vector_alloc(dim);
	design_mul0(1.0, TRANS_TRANS, design, isend, cm->p0, 0.0, cm->xbar0);

	return true;
}

static struct cohort_model *cohort_model_alloc(const struct design *design,
					       ssize_t isend,
					       const struct vector *coefs)
{
	struct cohort_model *cm;

	if ((cm = malloc(sizeof(*cm)))) {
		if (cohort_model_init(cm, design, isend, coefs))
			return cm;
		free(cm);
	}
	return NULL;
}

static void cohort_model_deinit(struct cohort_model *cm)
{
	assert(cm);
	vector_free(cm->xbar0);
	vector_free(cm->log_p0);
	vector_free(cm->p0);

}

static void cohort_model_free(struct cohort_model *cm)
{
	if (cm) {
		cohort_model_deinit(cm);
		free(cm);
	}
}

static void cohort_models_deinit(struct intmap *cohort_models)
{
	struct intmap_iter it;
	struct cohort_model *cm;

	intmap_iter_init(cohort_models, &it);
	while (intmap_iter_advance(cohort_models, &it)) {
		cm = *(struct cohort_model **)intmap_iter_current(cohort_models,
								  &it);
		cohort_model_deinit(cm);
	}
	intmap_iter_deinit(cohort_models, &it);
	intmap_deinit(cohort_models);
}

static bool insert_cohort_model(struct intmap *cohort_models,
				const struct design *design,
				ssize_t isend, const struct vector *coefs)
{
	const struct actors *senders = design_senders(design);
	intptr_t c = (intptr_t)actors_cohort(senders, isend);
	struct intmap_pos pos;
	struct cohort_model *cm;

	if (intmap_find(cohort_models, c, &pos))
		return true;

	if ((cm = cohort_model_alloc(design, isend, coefs))) {
		if (intmap_insert(cohort_models, &pos, &cm)) {
			return true;
		}

		cohort_model_free(cm);
	}
	return false;
}

static bool cohort_models_init(struct intmap *cohort_models,
			       const struct design *design,
			       const struct vector *coefs)
{
	ssize_t i, nsender = design_nsender(design);

	if (!intmap_init(cohort_models, sizeof(struct cohort_model *),
			 alignof(struct cohort_model *)))
		 return false;

	/* We have to loop over all senders, not over all cohorts, because we
	 * need a representative sender for each group when we call
	 * iproc_design_ctx_new(design, NULL, i).
	 */
	for (i = 0; i < nsender; i++) {
		if (!insert_cohort_model(cohort_models, design, i, coefs))
			goto fail;
	}
	return true;
fail:
	cohort_models_deinit(cohort_models);
	return false;
}

struct cohort_model *iproc_model_send_group(iproc_model * model, ssize_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	const struct design *design = iproc_model_design(model);
	const struct actors *senders = design_senders(design);
	intptr_t c = (intptr_t)actors_cohort(senders, isend);
	return *(struct cohort_model **)intmap_lookup(&model->cohort_models, c);
}

static void iproc_model_ctx_free_dealloc(iproc_model_ctx * ctx)
{
	if (ctx) {
		svector_free(ctx->dxbar);
		svector_free(ctx->dp);
		svector_free(ctx->deta);
		free(ctx);
	}
}

static void iproc_model_free(iproc_model * model)
{
	if (model) {
		ssize_t i, n = list_count(&model->ctxs);

		for (i = 0; i < n; i++) {
			iproc_model_ctx *ctx =
			    *(iproc_model_ctx **) list_item(&model->ctxs,
							    i);
			iproc_model_ctx_free_dealloc(ctx);
		}
		list_deinit(&model->ctxs);
		cohort_models_deinit(&model->cohort_models);
		vector_free(model->coefs);
		design_free(model->design);
		free(model);
	}
}

iproc_model *iproc_model_new(struct design *design,
			     struct vector *coefs, bool has_loops)
{
	assert(design);
	assert(coefs);
	assert(design_dim(design) == vector_dim(coefs));
	assert(design_nreceiver(design) > 0);
	assert(!has_loops || design_nreceiver(design) > 1);

	iproc_model *model = malloc(sizeof(*model));
	model->design = design_ref(design);
	model->coefs = vector_alloc_copy(coefs);
	model->has_loops = has_loops;
	cohort_models_init(&model->cohort_models, design, coefs);
	list_init(&model->ctxs, sizeof(iproc_model_ctx *));
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

struct design *iproc_model_design(iproc_model * model)
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

ssize_t iproc_model_nsender(iproc_model * model)
{
	assert(model);
	struct design *design = iproc_model_design(model);
	return design_nsender(design);
}

ssize_t iproc_model_nreceiver(iproc_model * model)
{
	assert(model);
	struct design *design = iproc_model_design(model);
	return design_nreceiver(design);
}

ssize_t iproc_model_dim(iproc_model * model)
{
	assert(model);
	struct design *design = iproc_model_design(model);
	return design_dim(design);
}

struct vector *iproc_model_logprobs0(iproc_model * model, ssize_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_group_model *group = iproc_model_send_group(model, isend);
	assert(group);
	return group->log_p0;
}

struct vector *iproc_model_probs0(iproc_model * model, ssize_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_group_model *group = iproc_model_send_group(model, isend);
	assert(group);
	return group->p0;
}

struct vector *iproc_model_mean0(iproc_model * model, ssize_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_group_model *group = iproc_model_send_group(model, isend);
	assert(group);
	return group->xbar0;
}
