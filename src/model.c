#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "logsumexp.h"
#include "model.h"
#include "util.h"

static void
compute_weight_changes(const struct frame *f, ssize_t isend,
		       const struct vector *coefs,
		       const struct vector *log_p0,
		       struct svector *deta, double *gamma, double *log_gamma)
{
	bool has_loops = design_loops(f->design);
	
	/* compute the changes in weights */
	frame_recv_dmul(1.0, TRANS_NOTRANS, f, isend, coefs, 0.0, deta);
	if (!has_loops) {
		svector_set_item(deta, isend, -INFINITY);
	}
	
	/* compute the scale for the weight differences */
	double lwmax = svector_max(deta);
	double logscale = MAX(0.0, lwmax);
	double invscale = exp(-logscale);
	
	struct logsumexp pos, neg;
	logsumexp_init(&pos);
	logsumexp_init(&neg);
	
	/* treat the initial sum of weights as the first positive difference */
	logsumexp_insert(&pos, -logscale);
	
	/* compute the log sums of the positive and negative differences in weights */
	struct svector_iter it;
	SVECTOR_FOREACH(it, deta) {
		//      for (i = 0; i < nnz; i++) {
		ssize_t jrecv = SVECTOR_IDX(it);
		double lp0 = *vector_item_ptr(log_p0, jrecv);
		double dlw = SVECTOR_VAL(it);
		double log_abs_dw;
		
		/* When w > w0:
		 *   w - w0      = w0 * [exp(dlw) - 1];
		 *               = w0 * {exp[dlw - log(scale)] - 1/scale} * scale
		 *
		 *   log(w - w0) = log(w0) + log{exp[(dlw) - log(scale)] - 1/scale} + log(scale)
		 */
		if (dlw >= 0) {
			log_abs_dw = lp0 + log(exp(dlw - logscale) - invscale);
			logsumexp_insert(&pos, log_abs_dw);
		} else {
			log_abs_dw = lp0 + log(invscale - exp(dlw - logscale));
			logsumexp_insert(&neg, log_abs_dw);
		}
	}
	
	double log_sum_abs_dw_p = logsumexp_value(&pos);
	double log_sum_abs_dw_n = logsumexp_value(&neg);
	
	if (log_sum_abs_dw_n > log_sum_abs_dw_p) {
		/* The sum of the weights is positive, so this only happens as a
		 * result of numerical errors */
		log_sum_abs_dw_n = log_sum_abs_dw_p;
		printf
		("\nWARNING: numerical errors in model-ctx.c (compute_new_logprobs)"
		 "\nPlease report a bug to the authors" "\n\n");
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
compute_active_probs(const struct vector *log_p0,
		     double log_gamma,
		     struct svector *deta, struct svector *p_active)
{
	ssize_t jrecv;
	double lp0_j, lp_j, deta_j, *p_j;
	struct svector_iter it;
	
	svector_assign_copy(p_active, deta);
	
	SVECTOR_FOREACH(it, p_active) {
		jrecv = SVECTOR_IDX(it);
		lp0_j = *vector_item_ptr(log_p0, jrecv);
		p_j = SVECTOR_PTR(it);
		deta_j = *p_j;
		lp_j = MIN(0.0, log_gamma + lp0_j + deta_j);
		*p_j = exp(lp_j);
	}
}

static void
compute_prob_diffs(const struct vector *p0, double gamma, struct svector *p_active)
{
	ssize_t jrecv;
	double p0_j, dp_j, *p_j;
	struct svector_iter it;
	
	SVECTOR_FOREACH(it, p_active) {
		jrecv = SVECTOR_IDX(it);
		p0_j = *vector_item_ptr(p0, jrecv);
		p_j = SVECTOR_PTR(it);
		dp_j = *p_j - gamma * p0_j;
		*p_j = dp_j;
		
	}
}


static void cohort_model_init(struct cohort_model *cm,
			      const struct design *design,
			      ssize_t isend, const struct vector *coefs)
{
	assert(cm);
	assert(design);
	assert(0 <= isend && isend < design_send_count(design));
	assert(coefs);
	
	ssize_t nreceiver = design_recv_count(design);
	ssize_t dim = design_recv_dim(design);

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
	design_recv_mul0(1.0, TRANS_NOTRANS, design, isend, coefs, 0.0, &cm->log_p0);
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
	design_recv_mul0(1.0, TRANS_TRANS, design, isend, &cm->p0, 0.0, &cm->xbar0);

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
	assert(0 <= isend && isend < model_sender_count(model));
	
	const struct design *design = model_design(model);
	const struct actors *senders = design_senders(design);
	intptr_t c = (intptr_t)actors_cohort(senders, isend);
	struct cohort_model *cm = intmap_item(&model->cohort_models, c);
	
	assert(cm);
	return cm;
}

static void send_model_clear(struct send_model *sm)
{
	assert(sm);
	
	const struct model *m = sm->model;
	ssize_t isend = sm->isend;
	const struct cohort_model *cm = model_cohort_model(m, isend);
	
	svector_clear(&sm->deta);	
	svector_clear(&sm->dp);
	svector_clear(&sm->dxbar);
	
	if (!design_loops(model_design(m))) {
		double p0 = vector_item(&cm->p0, isend);
		
		sm->gamma = 1.0 / (1.0 - p0);
		sm->log_gamma = -log1p(-p0);
		
		svector_set_item(&sm->deta, isend, -INFINITY);
		svector_set_item(&sm->dp, isend, -p0 / (1.0 - p0));
	} else {
		sm->gamma = 1.0;
		sm->log_gamma = 0.0;
	}
	sm->cached = true;
}

static void send_model_update(struct send_model *sm, const struct frame *f)
{
	assert(sm);
	assert(f);
	
	svector_clear(&sm->deta);	
	svector_clear(&sm->dp);
	svector_clear(&sm->dxbar);
	
	struct model *model = sm->model;
	ssize_t isend = sm->isend;
	const struct vector *coefs = model_coefs(model);
	const struct vector *log_p0 = model_logprobs0(model, isend);
	const struct vector *p0 = model_probs0(model, isend);
	
	compute_weight_changes(f, isend, coefs, log_p0,
			       &sm->deta, &sm->gamma, &sm->log_gamma);
	compute_active_probs(log_p0, sm->log_gamma, &sm->deta, &sm->dp);
	frame_recv_dmuls(1.0, TRANS_TRANS, f, isend, &sm->dp, 0.0, &sm->dxbar);
	compute_prob_diffs(p0, sm->gamma, &sm->dp);
	sm->cached = true;
}

static void send_model_init(struct send_model *sm, struct model *m, ssize_t isend)
{
	assert(sm);
	assert(m);
	assert(0 <= isend && isend < model_sender_count(m));
	
	sm->model = m;
	sm->isend = isend;
	sm->cohort = model_cohort_model(m, isend);
	svector_init(&sm->deta, model_receiver_count(m));
	svector_init(&sm->dp, model_receiver_count(m));
	svector_init(&sm->dxbar, model_dim(m));

	send_model_clear(sm);
}

static void send_model_deinit(struct send_model *sm)
{
	assert(sm);
	svector_deinit(&sm->dxbar);
	svector_deinit(&sm->dp);
	svector_deinit(&sm->deta);
}


static struct send_model *model_send_model_raw(struct model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < model_sender_count(m));
	
	struct intmap_pos pos;
	struct send_model *sm;
	
	if (!(sm = intmap_find(&m->send_models, isend, &pos))) {
		sm = intmap_insert(&m->send_models, &pos, NULL);
		send_model_init(sm, m, isend);
	}
	
	assert(sm);
	return sm;
}

void model_init(struct model *model, struct design *design,
		const struct vector *coefs)
{
	assert(model);
	assert(design);
	assert(coefs);
	assert(design_recv_dim(design) == vector_dim(coefs));
	assert(design_recv_count(design) > 0);
	assert(!design_loops(design) || design_recv_count(design) > 1);
	
	model->design = design_ref(design);
	vector_init_copy(&model->coefs, coefs);
	cohort_models_init(&model->cohort_models, design, coefs);
	intmap_init(&model->send_models, sizeof(struct send_model),
		    alignof(struct send_model));
	refcount_init(&model->refcount);
}

void model_deinit(struct model *model)
{
	assert(model);

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
	assert(design_recv_dim(design) == vector_dim(coefs));
	assert(design_recv_count(design) > 0);
	assert(!design_loops(design) || design_recv_count(design) > 1);

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
	return design_send_count(design);
}

ssize_t model_receiver_count(const struct model *model)
{
	assert(model);
	struct design *design = model_design(model);
	return design_recv_count(design);
}

ssize_t model_dim(const struct model *model)
{
	assert(model);
	struct design *design = model_design(model);
	return design_recv_dim(design);
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

void model_clear(struct model *m)
{
	assert(m);
	
	struct intmap_iter it;
	
	INTMAP_FOREACH(it, &m->send_models) {
		send_model_clear(INTMAP_VAL(it));
	}
}

static void process_recv_var_event(struct model *m,
				   const struct frame_event *e)
{
	assert(m);
	assert(e->type == RECV_VAR_EVENT);
	
	const struct recv_var_event_meta *meta = &e->meta.recv_var;
	ssize_t isend = meta->item.isend;
	struct send_model *sm = model_send_model_raw(m, isend);
	sm->cached = false;

	/*
	const struct cohort_model *cm = model_cohort_model(m, isend);
	ssize_t jrecv = meta->item.jrecv;
	struct svector_pos pos;
	double *ptr;
	
	ssize_t index = meta->index;
	double beta = vector_item(model_coefs(m), index);
	double dx = meta->delta;
	double deta_j = beta * dx;
	double p0 = sm->gamma * vector_item(&cm->p0, jrecv);
	
	ptr = svector_find(&sm->deta, jrecv, &pos);
	if (!ptr) {
		svector_insert(&sm->deta, &pos, deta_j);
	} else {
		p0 += svector_item(&sm->dp, jrecv);
		*ptr += deta_j;
		if (*ptr == 0.0) {
			svector_remove_at(&sm->deta, &pos);
		}
	}*/
}

static void model_update_with(struct model *m, const struct frame_event *e)
{
	assert(m);
	assert(e);
	
	switch (e->type) {
	case RECV_VAR_EVENT:
		process_recv_var_event(m, e);
		break;
	default:
		break; /* pass */
	}
}

void model_update(struct model *m, const struct frame *f)
{
	assert(m);
	assert(f);
	
	const struct frame_event *e;
	
	ARRAY_FOREACH(e, &f->events) {
		if (e->type & (SEND_VAR_EVENT | RECV_VAR_EVENT)) {
			model_update_with(m, e);
		}
	}
}

struct send_model *model_send_model(struct model *m, const struct frame *f, ssize_t isend)
{
	assert(m);
	assert(f);
	
	struct send_model *sm = model_send_model_raw(m, isend);
	if (!sm->cached) {
		send_model_update(sm, f);
	}
	return sm;
}

ssize_t send_model_receiver_count(const struct send_model *sm)
{
	assert(sm);
	struct model *model = sm->model;
	ssize_t nreceiver = model_receiver_count(model);
	return nreceiver;
}

/*
 *         log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double send_model_logprob(const struct send_model *sm, ssize_t jrecv)
{
	assert(sm);
	assert(0 <= jrecv && jrecv < send_model_receiver_count(sm));
	
	/*
	 double gamma = ctx->gamma;
	 double p0 = vector_item(ctx->group->p0, jrecv);
	 double dp = svector_item(ctx->dp, jrecv);
	 double p = gamma * p0 + dp;
	 p = MAX(0.0, MIN(1.0, p));
	 return log(p);
	 */
	
	double log_gamma = sm->log_gamma;
	double log_p0 = *vector_item_ptr(&sm->cohort->log_p0, jrecv);
	double deta = svector_item(&sm->deta, jrecv);
	double log_p = log_gamma + log_p0 + deta;
	return MIN(log_p, 0.0);
}

double send_model_prob(const struct send_model *sm, ssize_t jrecv)
{
	assert(sm);
	assert(0 <= jrecv && jrecv < send_model_receiver_count(sm));
	
	double lp = send_model_logprob(sm, jrecv);
	double p = exp(lp);
	return p;
}

void send_model_get_logprobs(const struct send_model *sm, struct vector *logprobs)
{
	assert(sm);
	assert(logprobs);
	assert(vector_dim(logprobs) == send_model_receiver_count(sm));
	
	ssize_t j, n = send_model_receiver_count(sm);
	
	for (j = 0; j < n; j++) {
		double lp = send_model_logprob(sm, j);
		*vector_item_ptr(logprobs, j) = lp;
	}
}

void send_model_get_probs(const struct send_model *sm, struct vector *probs)
{
	assert(sm);
	assert(probs);
	assert(vector_dim(probs) == send_model_receiver_count(sm));
	
	send_model_get_logprobs(sm, probs);
	vector_exp(probs);
}
