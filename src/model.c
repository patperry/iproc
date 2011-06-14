#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "model.h"
#include "util.h"

/* DEPRECATED */
static void iproc_model_ctx_free_dealloc(iproc_model_ctx * ctx)
{
	if (ctx) {
		svector_deinit(&ctx->dxbar);
		svector_deinit(&ctx->dp);
		svector_deinit(&ctx->deta);
		free(ctx);
	}
}




static void send_model_init(struct send_model *sm, struct model *m, ssize_t isend);
static void send_model_deinit(struct send_model *sm);

static void cohort_model_init(struct cohort_model *cm,
			      const struct design *design,
			      ssize_t isend, const struct vector *coefs)
{
	assert(cm);
	assert(design);
	assert(0 <= isend && isend < design_nsender(design));
	assert(coefs);
	
	ssize_t nreceiver = design_nreceiver(design);
	ssize_t dim = design_dim(design);

	/* The probabilities are p[i] = w[i] / sum(w[j]), so
	 * log(p[i]) = log(w[i]) - log(sum(w[j])).
	 */

	/* NOTE: should we worry about overflow in the multiplication? 
	 * It is possible to gaurd against said overflow by computing a
	 * scaled version of eta0: scale the coefficient vector before the
	 * multiplication, then unscale when computing p0.  This
	 * shouldn't be necessary in most (all?) real-world situations.
	 */
	
	/* log_p0, eta0 */
	vector_init(&cm->log_p0, nreceiver);
	design_mul0(1.0, TRANS_NOTRANS, design, isend, coefs, 0.0, &cm->log_p0);
#ifndef NDEGUG
	vector_init_copy(&cm->eta0, &cm->log_p0);
#endif
	double max = vector_max(&cm->log_p0); /* protect against overflow */
	vector_shift(&cm->log_p0, -max);
	double lse = vector_log_sum_exp(&cm->log_p0);
	vector_shift(&cm->log_p0, -lse);
	
	/* log_W0 */
	cm->log_W0 = lse + max;
	
	/* p0 */
	vector_init_copy(&cm->p0, &cm->log_p0);
	vector_exp(&cm->p0);
	
	/* xbar0 */
	vector_init(&cm->xbar0, dim);
	design_mul0(1.0, TRANS_TRANS, design, isend, &cm->p0, 0.0, &cm->xbar0);

#ifndef NDEBUG
	/* W0 */
	cm->W0 = exp(cm->log_W0);

	/* w0 */
	vector_init_copy(&cm->w0, &cm->p0);
	vector_scale(&cm->w0, cm->W0);
#endif

}

static void cohort_model_deinit(struct cohort_model *cm)
{
	assert(cm);
#ifndef NDEBUG
	vector_deinit(&cm->w0);
	vector_deinit(&cm->eta0);	
#endif
	vector_deinit(&cm->xbar0);
	vector_deinit(&cm->log_p0);
	vector_deinit(&cm->p0);
}

static void cohort_models_init(struct intmap *cohort_models,
			       const struct design *design,
			       const struct vector *coefs)
{
	const struct actors *senders = design_senders(design);
	ssize_t isend, nsend = actors_count(senders);
	struct intmap_pos pos;
	struct cohort_model *cm;
	intptr_t c;
	
	intmap_init(cohort_models, sizeof(struct cohort_model),
		    alignof(struct cohort_model));

	/* We have to loop over all senders, not over all cohorts, because we
	 * need a representative sender for each group when we call
	 * iproc_design_ctx_new(design, NULL, i).
	 */
	for (isend = 0; isend < nsend; isend++) {
		c = (intptr_t)actors_cohort(senders, isend);
		if (intmap_find(cohort_models, c, &pos))
			continue;
		
		cm = intmap_insert(cohort_models, &pos, NULL);
		cohort_model_init(cm, design, isend, coefs);
	}
}

static void cohort_models_deinit(struct intmap *cohort_models)
{
	struct intmap_iter it;
	struct cohort_model *cm;
	
	INTMAP_FOREACH(it, cohort_models) {
		cm = INTMAP_VAL(it);
		cohort_model_deinit(cm);
	}
	intmap_deinit(cohort_models);
}

static struct cohort_model *model_cohort_model(const struct model *model,
					       ssize_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < model_sender_count(model));
	
	const struct design *design = model_design(model);
	const struct actors *senders = design_senders(design);
	intptr_t c = (intptr_t)actors_cohort(senders, isend);
	struct cohort_model *cm = intmap_item(&model->cohort_models, c);
	
	assert(cm);
	return cm;
}

void model_init(struct model *model, struct design *design,
		const struct vector *coefs)
{
	assert(model);
	assert(design);
	assert(coefs);
	assert(design_dim(design) == vector_dim(coefs));
	assert(design_nreceiver(design) > 0);
	assert(!design_loops(design) || design_nreceiver(design) > 1);
	
	model->design = design_ref(design);
	vector_init_copy(&model->coefs, coefs);
	cohort_models_init(&model->cohort_models, design, coefs);
	intmap_init(&model->send_models, sizeof(struct send_model),
		    alignof(struct send_model));

	
	array_init(&model->ctxs, sizeof(iproc_model_ctx *));
	refcount_init(&model->refcount);
}

void model_deinit(struct model *model)
{
	assert(model);

	/* DEPRECATED */
	iproc_model_ctx *ctx;
	ARRAY_FOREACH(ctx, &model->ctxs) {
		iproc_model_ctx_free_dealloc(ctx);
	}
	array_deinit(&model->ctxs);
	
	refcount_deinit(&model->refcount);
	
	struct intmap_iter it;
	INTMAP_FOREACH(it, &model->send_models) {
		send_model_deinit(INTMAP_VAL(it));
	}
	
	intmap_deinit(&model->send_models);
	cohort_models_deinit(&model->cohort_models);
	vector_deinit(&model->coefs);
	design_free(model->design);
}



struct model *model_alloc(struct design *design,
			  const struct vector *coefs)
{
	assert(design);
	assert(coefs);
	assert(design_dim(design) == vector_dim(coefs));
	assert(design_nreceiver(design) > 0);
	assert(!design_loops(design) || design_nreceiver(design) > 1);

	struct model *model = xcalloc(1, sizeof(*model));
	model_init(model, design, coefs);
	return model;
}


struct model *model_ref(struct model *model)
{
	assert(model);
	refcount_get(&model->refcount);
	return model;
}

void model_free(struct model *model)
{
	if (model && refcount_put(&model->refcount, NULL)) {
		refcount_get(&model->refcount);
		model_deinit(model);
		xfree(model);
	}
}



void send_model_init(struct send_model *sm, struct model *m, ssize_t isend)
{
	assert(sm);
	assert(m);
	assert(0 <= isend && isend < model_sender_count(m));
	
	sm->model = m;
	sm->isend = isend;
	svector_init(&sm->deta, model_receiver_count(m));
	svector_init(&sm->dp, model_receiver_count(m));
	svector_init(&sm->dxbar, model_receiver_count(m));
}

void send_model_deinit(struct send_model *sm)
{
	assert(sm);
	svector_deinit(&sm->dxbar);
	svector_deinit(&sm->dp);
	svector_deinit(&sm->deta);
}


struct design *model_design(const struct model *model)
{
	assert(model);
	return model->design;
}

struct vector *model_coefs(const struct model *model)
{
	assert(model);
	return &((struct model *)model)->coefs;
}

ssize_t model_sender_count(const struct model *model)
{
	assert(model);
	struct design *design = model_design(model);
	return design_nsender(design);
}

ssize_t model_receiver_count(const struct model *model)
{
	assert(model);
	struct design *design = model_design(model);
	return design_nreceiver(design);
}

ssize_t model_dim(const struct model *model)
{
	assert(model);
	struct design *design = model_design(model);
	return design_dim(design);
}

struct vector *model_logprobs0(const struct model *model, ssize_t isend)
{
	assert(model);
	assert(0 <= isend && isend < model_sender_count(model));

	struct cohort_model *cm = model_cohort_model(model, isend);
	return &cm->log_p0;
}

struct vector *model_probs0(const struct model *model, ssize_t isend)
{
	assert(model);
	assert(0 <= isend && isend < model_sender_count(model));

	struct cohort_model *cm = model_cohort_model(model, isend);
	return &cm->p0;
}

struct vector *model_mean0(const struct model *model, ssize_t isend)
{
	assert(model);
	assert(0 <= isend && isend < model_sender_count(model));

	struct cohort_model *cm = model_cohort_model(model, isend);
	return &cm->xbar0;
}
