#include "port.h"

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "model.h"
#include "util.h"

static void
compute_logprobs0(const iproc_design * design,
		  int64_t isend,
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

static bool cohort_model_init (struct cohort_model *cm,
			       const struct design *design,
			       ssize_t isend,
			       const struct vector *coefs)
{
	ssize_t nreceiver = iproc_design_nreceiver(design);
	ssize_t dim = iproc_design_dim(design);

	/* compute initial log(probs) */
	cm->log_p0 = vector_new(nreceiver);
	compute_logprobs0(design, isend, coefs, cm->log_p0,
			  &cm->log_W0);
	
	/* compute initial probs */
	cm->p0 = vector_new_copy(cm->log_p0);
	vector_exp(cm->p0);
	
	/* compute initial covariate mean */
	cm->xbar0 = vector_new(dim);
	iproc_design_mul0(1.0, IPROC_TRANS_TRANS, design, isend, cm->p0,
			  0.0, cm->xbar0);
	
	return true;
}

static void cohort_model_deinit (struct cohort_model *cm)
{
	assert(cm);
	vector_free(cm->xbar0);
	vector_free(cm->log_p0);
	vector_free(cm->p0);
	
}


static bool cohort_models_init(struct intmap *cohort_models,
			       const iproc_design * design,
			       const struct vector *coefs)
{
	ssize_t i, nsender = iproc_design_nsender(design);
	const iproc_actors *senders = iproc_design_senders(design);
	struct intmap_pos pos;
	struct cohort_model *cm;
	intptr_t c;

	if (!intmap_init(cohort_models, struct cohort_model))
		return false;

	/* We have to loop over all senders, not over all cohorts, because we
	 * need a representative sender for each group when we call
	 * iproc_design_ctx_new(design, NULL, i).
	 */
	for (i = 0; i < nsender; i++) {
		c = iproc_actors_group(senders, i);
		if (!intmap_find(cohort_models, c, &pos)) {
			cm = intmap_insert(cohort_models, &pos, NULL);
			cohort_model_init(cm, design, i, coefs);
		}
	}
	return true;
}

static void cohort_models_deinit(struct intmap *cohort_models)
{
	struct cohort_model *cm = intmap_vals_begin(cohort_models);
	struct cohort_model *end = intmap_vals_end(cohort_models);
	
	for (; cm < end; cm++) {
		cohort_model_deinit(cm);
	}

	intmap_deinit(cohort_models);
}

struct cohort_model *iproc_model_send_group(iproc_model * model, int64_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	const iproc_design *design = iproc_model_design(model);
	const iproc_actors *senders = iproc_design_senders(design);
	int64_t c = iproc_actors_group(senders, isend);
	return intmap_lookup(&model->cohort_models, c);
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
		cohort_models_deinit(&model->cohort_models);
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
	cohort_models_init(&model->cohort_models, design, coefs);
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
