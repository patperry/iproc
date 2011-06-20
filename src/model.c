#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "logsumexp.h"
#include "model.h"
#include "util.h"

DEFINE_COMPARE_FN(ssize_compare, ssize_t)

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
compute_prob_diffs(const struct vector *p0, double gamma,
		   struct svector *p_active)
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
	design_recv_mul0(1.0, TRANS_NOTRANS, design, isend, coefs, 0.0,
			 &cm->log_p0);
#ifndef NDEGUG
	vector_init_copy(&cm->eta0, &cm->log_p0);
#endif
	double max = vector_max(&cm->log_p0);	/* protect against overflow */
	vector_shift(&cm->log_p0, -max);
	double lse = vector_log_sum_exp(&cm->log_p0);
	vector_shift(&cm->log_p0, -lse);

	/* log_W0 */
	cm->log_W0 = lse + max;

	/* p0 */
	vector_init_copy(&cm->p0, &cm->log_p0);
	vector_exp(&cm->p0);

	/* mean0 */
	vector_init(&cm->mean0, dim);
	design_recv_mul0(1.0, TRANS_TRANS, design, isend, &cm->p0, 0.0,
			 &cm->mean0);
	
	/* imat0 */
	struct vector y;
	struct svector ej;
	double pj;
	ssize_t jrecv;	
	
	vector_init(&y, dim);
	svector_init(&ej, nreceiver);
	matrix_init(&cm->imat0, dim, dim);

	for (jrecv = 0; jrecv < nreceiver; jrecv++) {
		vector_assign_copy(&y, &cm->mean0);
		svector_set_basis(&ej, jrecv);
		design_recv_muls0(1.0, TRANS_TRANS, design, isend, &ej, -1.0, &y);
		pj = vector_item(&cm->p0, jrecv);
		
		matrix_update1(&cm->imat0, pj, &y, &y);
	}
	svector_deinit(&ej);
	vector_deinit(&y);

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
	matrix_deinit(&cm->imat0);	
	vector_deinit(&cm->mean0);
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

static void recv_model_clear(struct recv_model *rm)
{
	assert(rm);

	const struct model *m = rm->model;
	ssize_t isend = rm->isend;
	const struct cohort_model *cm = model_cohort_model(m, isend);

	svector_clear(&rm->deta);
	array_clear(&rm->active);
	vector_deinit(&rm->dp);
	vector_fill(&rm->mean_dx, 0.0);
	
	matrix_deinit(&rm->dp2);
	matrix_fill(&rm->var_dx, 0.0);	

	if (!design_loops(model_design(m))) {
		double p0 = vector_item(&cm->p0, isend);
		double dp = -p0 / (1.0 - p0);

		rm->gamma = 1.0 / (1.0 - p0);
		rm->log_gamma = -log1p(-p0);

		svector_set_item(&rm->deta, isend, -INFINITY);
		
		array_add(&rm->active, &isend);
		vector_init(&rm->dp, 1);
		vector_set_item(&rm->dp, 0, dp);
		matrix_init(&rm->dp2, 1, 1);
		matrix_set_item(&rm->dp2, 0, 0, dp * (1 - dp));
	} else {
		rm->gamma = 1.0;
		rm->log_gamma = 0.0;
		vector_init(&rm->dp, 0);
		matrix_init(&rm->dp2, 0, 0);
	}
	rm->cached = true;
}

static void recv_model_update(struct recv_model *rm, const struct frame *f)
{
	assert(rm);
	assert(f);

	//vector_fill(&rm->dp, 0.0);
	vector_fill(&rm->mean_dx, 0.0);

	ssize_t isend = rm->isend;
	//const struct vector *log_p0 = recv_model_logprobs0(rm);
	const struct vector *probs0 = recv_model_probs0(rm);
	double gamma = rm->gamma;
	
	//compute_weight_changes(f, rm->isend, model_coefs(rm->model), log_p0,
	//		       &rm->deta, &rm->gamma, &rm->log_gamma);

	
	ssize_t i, n = array_count(&rm->active);
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&rm->active, i);

		/* mean */		
		double p = recv_model_prob(rm, jrecv);
		double p0 = vector_item(probs0, jrecv);
		double dp = p - gamma * p0;
		vector_set_item(&rm->dp, i, dp);
		
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		vector_axpy(p, dx, &rm->mean_dx);
	}

	/* var */
	matrix_fill(&rm->dp2, 0.0);	
	matrix_fill(&rm->var_dx, 0.0);
	struct vector y;
	double ptot = 0.0;

	matrix_update1(&rm->dp2, -1.0, &rm->dp, &rm->dp);
	
	vector_init(&y, vector_dim(&rm->mean_dx));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&rm->active, i);
		double dp = vector_item(&rm->dp, i);
		
		*matrix_item_ptr(&rm->dp2, i, i) += dp;
		
		double p = recv_model_prob(rm, jrecv);
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		
		vector_assign_copy(&y, dx);
		vector_sub(&y, &rm->mean_dx);
		
		matrix_update1(&rm->var_dx, p, &y, &y);
		ptot += p;
	}
	matrix_update1(&rm->var_dx, 1 - ptot, &rm->mean_dx, &rm->mean_dx);
	vector_deinit(&y);

	

	//compute_active_probs(log_p0, rm->log_gamma, &rm->deta, &rm->dp);
	//frame_recv_dmuls(1.0, TRANS_TRANS, f, isend, &rm->dp, 0.0, &rm->mean_dx);
	//compute_prob_diffs(p0, rm->gamma, &rm->dp);
	rm->cached = true;	
}

static void recv_model_init(struct recv_model *rm, struct model *m,
			    ssize_t isend)
{
	assert(rm);
	assert(m);
	assert(0 <= isend && isend < model_sender_count(m));

	ssize_t dim_dx = design_recv_dyn_dim(model_design(m));
	
	rm->model = m;
	rm->isend = isend;
	rm->cohort = model_cohort_model(m, isend);
	svector_init(&rm->deta, model_receiver_count(m));
	array_init(&rm->active, sizeof(ssize_t));
	vector_init(&rm->dp, 0);
	vector_init(&rm->mean_dx, dim_dx);
	
	matrix_init(&rm->dp2, 0, 0);
	matrix_init(&rm->var_dx, dim_dx, dim_dx);

	recv_model_clear(rm);
}

static void recv_model_deinit(struct recv_model *rm)
{
	assert(rm);
	matrix_deinit(&rm->var_dx);
	matrix_deinit(&rm->dp2);
	vector_deinit(&rm->mean_dx);
	vector_deinit(&rm->dp);
	array_deinit(&rm->active);
	svector_deinit(&rm->deta);
}

static struct recv_model *model_recv_model_raw(struct model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < model_sender_count(m));

	struct intmap_pos pos;
	struct recv_model *rm;

	if (!(rm = intmap_find(&m->recv_models, isend, &pos))) {
		rm = intmap_insert(&m->recv_models, &pos, NULL);
		recv_model_init(rm, m, isend);
	}

	assert(rm);
	return rm;
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
	intmap_init(&model->recv_models, sizeof(struct recv_model),
		    alignof(struct recv_model));
	refcount_init(&model->refcount);
}

void model_deinit(struct model *model)
{
	assert(model);

	refcount_deinit(&model->refcount);

	struct intmap_iter it;
	INTMAP_FOREACH(it, &model->recv_models) {
		recv_model_deinit(INTMAP_VAL(it));
	}

	intmap_deinit(&model->recv_models);
	cohort_models_deinit(&model->cohort_models);
	vector_deinit(&model->coefs);
	design_free(model->design);
}

struct model *model_alloc(struct design *design, const struct vector *coefs)
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

struct vector *recv_model_logprobs0(const struct recv_model *model)
{
	assert(model);
	return &model->cohort->log_p0;
}

struct vector *recv_model_probs0(const struct recv_model *model)
{
	assert(model);
	return &model->cohort->p0;
}

struct vector *recv_model_mean0(const struct recv_model *model)
{
	assert(model);
	return &model->cohort->mean0;
}

struct matrix *recv_model_imat0(const struct recv_model *model)
{
	assert(model);
	return &model->cohort->imat0;
}

void model_clear(struct model *m)
{
	assert(m);

	struct intmap_iter it;

	INTMAP_FOREACH(it, &m->recv_models) {
		recv_model_clear(INTMAP_VAL(it));
	}
}


static void process_recv_var_event(struct model *m, const struct frame_event *e)
{
	assert(m);
	assert(e->type == RECV_VAR_EVENT);

	const struct recv_var_event_meta *meta = &e->meta.recv_var;
	struct recv_model *rm = model_recv_model_raw(m, meta->item.isend);
	const struct cohort_model *cm = rm->cohort;	
	const struct vector *coefs = model_coefs(m);	

	ssize_t jrecv = meta->item.jrecv;
	double *pdeta = svector_item_ptr(&rm->deta, jrecv);

	ssize_t index = meta->index;
	double dx = meta->delta;
	double deta = vector_item(coefs, index) * dx;

	double gamma0 = rm->gamma;
	double log_gamma0 = rm->log_gamma;
	// double p0 = gamma0 * vector_item(&cm->p0, jrecv) * exp(*pdeta);
	double log_p0 = log_gamma0 + vector_item(&cm->log_p0, jrecv) + (*pdeta);
	double p0 = exp(log_p0);
	double p = exp(log_p0 + deta);
	double dp = p - p0;
	double gamma = gamma0 / (1.0 + dp);
	double log_gamma = log_gamma0 - log1p(dp);
	
	rm->gamma = gamma;
	rm->log_gamma = log_gamma;
	*pdeta += deta;
	
	ssize_t pos = array_binary_search(&rm->active, &jrecv, ssize_compare);
	if (pos < 0) {
		pos = ~pos;
		array_insert(&rm->active, pos, &jrecv);
		
		ssize_t dim = vector_dim(&rm->dp) + 1;
		struct vector dp, dp0, dp1, dp0_old, dp1_old;
		
		vector_init(&dp, dim);
		dp0 = vector_slice(&dp, 0, pos);
		dp0_old = vector_slice(&rm->dp, 0, pos);
		vector_assign_copy(&dp0, &dp0_old);
		
		dp1 = vector_slice(&dp, pos + 1, dim - pos - 1);
		dp1_old = vector_slice(&rm->dp, pos, dim - pos - 1);
		vector_assign_copy(&dp1, &dp1_old);
		
		vector_deinit(&rm->dp);
		rm->dp = dp;
		
		struct matrix dp2, dst, src;
		matrix_init(&dp2, dim, dim);
		dst = matrix_slice(&dp2, 0, 0, pos, pos);
		src = matrix_slice(&rm->dp2, 0, 0, pos, pos);
		matrix_assign_copy(&dst, &src);
		
		dst = matrix_slice(&dp2, pos + 1, 0, dim - pos - 1, pos);
		src = matrix_slice(&rm->dp2, pos, 0, dim - pos - 1, pos);
		matrix_assign_copy(&dst, &src);

		dst = matrix_slice(&dp2, 0, pos + 1, pos, dim - pos - 1);
		src = matrix_slice(&rm->dp2, 0, pos, pos, dim - pos - 1);
		matrix_assign_copy(&dst, &src);

		dst = matrix_slice(&dp2, pos + 1, pos + 1, dim - pos - 1, dim - pos - 1);
		src = matrix_slice(&rm->dp2, pos, pos, dim - pos - 1, dim - pos - 1);
		matrix_assign_copy(&dst, &src);
		matrix_deinit(&rm->dp2);
		rm->dp2 = dp2;
	}
	
	rm->cached = false;	
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
		break;		/* pass */
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

struct recv_model *model_recv_model(const struct model *m, const struct frame *f,
				    ssize_t isend)
{
	assert(m);

	struct recv_model *rm = model_recv_model_raw((struct model *)m, isend);
	if (!rm->cached) {
		assert(f);
		recv_model_update(rm, f);
	}
	return rm;
}

ssize_t recv_model_count(const struct recv_model *rm)
{
	assert(rm);
	const struct model *model = rm->model;
	const struct design *design = model_design(model);
	return design_recv_count(design);
}

ssize_t recv_model_dim(const struct recv_model *rm)
{
	assert(rm);
	const struct model *model = rm->model;
	const struct design *design = model_design(model);
	return design_recv_dim(design);
}


/*
 *         log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double recv_model_logprob(const struct recv_model *rm, ssize_t jrecv)
{
	assert(rm);
	assert(0 <= jrecv && jrecv < recv_model_count(rm));

	/*
	   double gamma = ctx->gamma;
	   double p0 = vector_item(ctx->group->p0, jrecv);
	   double dp = svector_item(ctx->dp, jrecv);
	   double p = gamma * p0 + dp;
	   p = MAX(0.0, MIN(1.0, p));
	   return log(p);
	 */

	double log_gamma = rm->log_gamma;
	double log_p0 = *vector_item_ptr(&rm->cohort->log_p0, jrecv);
	double deta = svector_item(&rm->deta, jrecv);
	double log_p = log_gamma + log_p0 + deta;
	return MIN(log_p, 0.0);
}

double recv_model_prob(const struct recv_model *rm, ssize_t jrecv)
{
	assert(rm);
	assert(0 <= jrecv && jrecv < recv_model_count(rm));

	double lp = recv_model_logprob(rm, jrecv);
	double p = exp(lp);
	return p;
}

void recv_model_axpy_probs(double alpha,
			   const struct recv_model *rm,
			   struct vector *y)
{
	assert(rm);
	assert(y);
	assert(vector_dim(y) == recv_model_count(rm));

	ssize_t j, n = recv_model_count(rm);

	for (j = 0; j < n; j++) {
		double p = recv_model_prob(rm, j);
		*vector_item_ptr(y, j) += alpha * p;
	}
}

void recv_model_axpy_mean(double alpha,
			  const struct recv_model *rm,
			  struct vector *y)
{
	assert(rm);
	assert(y);
	assert(vector_dim(y) == recv_model_dim(rm));
	assert(rm->cached);
	
	const struct vector *mean0 = recv_model_mean0(rm);
	double gamma = rm->gamma;
	vector_axpy(alpha * gamma, mean0, y);
	
	const struct design *design = model_design(rm->model);
	ssize_t isend = rm->isend;
	
	ssize_t i, n = array_count(&rm->active);
	struct svector dp;
	svector_init(&dp, recv_model_count(rm));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&rm->active, i);
		double val = vector_item(&rm->dp, i);
		svector_set_item(&dp, jrecv, val);
	}
	
	design_recv_muls0(alpha, TRANS_TRANS, design, isend, &dp, 1.0, y);
	svector_deinit(&dp);
	
	const struct vector *mean_dx = &rm->mean_dx;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);
	vector_axpy(alpha, mean_dx, &ysub);
}

void recv_model_axpy_var(double alpha, const struct recv_model *rm, const struct frame *f, struct matrix *y)
{
	assert(rm);
	assert(y);
	assert(matrix_ncol(y) == recv_model_dim(rm));
	assert(matrix_nrow(y) == recv_model_dim(rm));	
	assert(rm->cached);
	
	const ssize_t dim = recv_model_dim(rm);
	const struct design *design = model_design(rm->model);
	const ssize_t isend = rm->isend;
	const struct vector *mean0 = recv_model_mean0(rm);
	const struct matrix *var0 = &rm->cohort->imat0;
	const double gamma = rm->gamma;
	
	const ssize_t n = array_count(&rm->active);
	ssize_t i, j, k;
	struct svector dp;
	svector_init(&dp, recv_model_count(rm));
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&rm->active, i);
		double val = vector_item(&rm->dp, i);
		svector_set_item(&dp, jrecv, val);
	}
	
	struct vector x0_dp;
	vector_init(&x0_dp, dim);
	design_recv_muls0(1.0, TRANS_TRANS, design, isend, &dp, 0.0, &x0_dp);

	struct matrix x0_dp2;
	struct svector dp2_j;
	
	matrix_init(&x0_dp2, dim, n);
	svector_init(&dp2_j, recv_model_count(rm));
	for (j = 0; j < n; j++) {
		svector_clear(&dp2_j);
		for (i = 0; i < n; i++) {
			ssize_t jrecv = *(ssize_t *)array_item(&rm->active, i);
			double val = matrix_item(&rm->dp2, i, j);
			svector_set_item(&dp2_j, jrecv, val);
		}
		struct vector dst = matrix_col(&x0_dp2, j);
		design_recv_muls0(1.0, TRANS_TRANS, design, isend, &dp2_j, 0.0, &dst);
	}
	
	matrix_axpy(alpha * gamma, var0, y);
	matrix_update1(y, alpha * gamma * (1 - gamma), mean0, mean0);
	matrix_update1(y, -alpha * gamma, mean0, &x0_dp);
	matrix_update1(y, -alpha * gamma, &x0_dp, mean0);
	
	struct svector x0_dp2_k;
	svector_init(&x0_dp2_k, recv_model_count(rm));
	for (k = 0; k < dim; k++) {
		svector_clear(&x0_dp2_k);
		for (j = 0; j < n; j++) {
			ssize_t jrecv = *(ssize_t *)array_item(&rm->active, j);
			double val = matrix_item(&x0_dp2, k, j);
			svector_set_item(&x0_dp2_k, jrecv, val);
		}
		struct vector dst = matrix_col(y, k);
		design_recv_muls0(alpha, TRANS_TRANS, design, isend, &x0_dp2_k, 1.0, &dst);
	}
	
	ssize_t dyn_off = design_recv_dyn_index(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);	
	struct matrix y_1 = matrix_slice(y, 0, dyn_off, dim, dyn_dim);
	struct matrix y1_ = matrix_slice(y, dyn_off, 0, dyn_dim, dim);
	struct matrix y11 = matrix_slice(y, dyn_off, dyn_off, dyn_dim, dyn_dim);

	matrix_axpy(alpha, &rm->var_dx, &y11);
	
	struct vector x0_j;	
	struct svector e_j;

	vector_init(&x0_j, dim);
	svector_init(&e_j, design_recv_count(design));
	
	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&rm->active, i);
		double p = recv_model_prob(rm, jrecv);
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);		
		
		svector_set_basis(&e_j, jrecv);
		design_recv_muls0(alpha, TRANS_TRANS, design, isend, &e_j, 0.0, &x0_j);
		
		matrix_update1(&y_1, alpha * p, &x0_j, dx);
		matrix_update1(&y1_, alpha * p, dx, &x0_j);
	}
	
	matrix_update1(&y_1, -alpha * gamma, mean0, &rm->mean_dx);
	matrix_update1(&y1_, -alpha * gamma, &rm->mean_dx, mean0);

	matrix_update1(&y_1, -alpha, &x0_dp, &rm->mean_dx);
	matrix_update1(&y1_, -alpha, &rm->mean_dx, &x0_dp);

	vector_deinit(&x0_j);
	svector_deinit(&e_j);
	svector_deinit(&x0_dp2_k);
	svector_deinit(&dp2_j);
	matrix_deinit(&x0_dp2);
	vector_deinit(&x0_dp);
	svector_deinit(&dp);
}

