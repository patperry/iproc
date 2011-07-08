#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "logsumexp.h"
#include "model.h"
#include "util.h"

DEFINE_COMPARE_FN(ssize_compare, ssize_t)

static void
compute_weight_changes(const struct frame *f, ssize_t isend,
		       const struct vector *recv_coefs,
		       const struct vector *eta0,
		       double max_eta0, double log_W0,
		       const struct svector *deta, double *scale,
		       double *gamma, double *log_W)
{
	//if (isend == 119) {
	//	printf("COMPUTE\n");
	//}
	
	const struct design *design = f->design;
	
	ssize_t jrecv, nrecv = design_recv_count(design);				
	struct svector_iter it;
	bool shrink = false;
	
	/* compute the maximum eta value */
	*scale = max_eta0;

	SVECTOR_FOREACH(it, deta) {
		ssize_t jrecv = SVECTOR_IDX(it);
		double eta0_j = vector_item(eta0, jrecv);
		double deta_j = SVECTOR_VAL(it);
		double eta_j = eta0_j + deta_j;
		
		if (eta0_j == max_eta0 && deta_j < 0)
			shrink = true;

		*scale = MAX(*scale, eta_j);
	}
	
fast:
	{
		double W = exp(log_W0 + (max_eta0 - *scale));
		bool found_max = false;
		
		SVECTOR_FOREACH(it, deta) {
			ssize_t jrecv = SVECTOR_IDX(it);
			double eta0_j = vector_item(eta0, jrecv);
			double deta_j = SVECTOR_VAL(it);
			double eta_j = eta0_j + deta_j;
			
			if (!found_max && eta_j == *scale) {
				found_max = true;
				W += -exp(eta0_j - *scale);
			} else {
				W += (exp(eta_j - *scale)
				      - exp(eta0_j - *scale));
			}
		}
			
		if (found_max) {
			*log_W = log1p(W);
		} else {
			*log_W = log(W);
		}
		
		if (*log_W > 0)
			goto out;
	}

accurate:
	{
		//fprintf(stderr, "!"); fflush(stderr);
		/* compute the new max_eta */
		if (shrink && *scale == max_eta0) {
			double max = -INFINITY;
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double eta0_j = vector_item(eta0, jrecv);
				double deta_j = svector_item(deta, jrecv);
				double eta_j = eta0_j + deta_j;
				
				max = MAX(max, eta_j);
			}
			*scale = max;
		}

		/* compute the new log_W */
		struct logsumexp lse;
		logsumexp_init(&lse);		

		for (jrecv = 0; jrecv < nrecv; jrecv++) {
			double eta0_j = vector_item(eta0, jrecv);
			double deta_j = svector_item(deta, jrecv);
			double eta_j = eta0_j + deta_j;
		
			logsumexp_insert(&lse, eta_j - *scale);
		}
		*log_W = logsumexp_value(&lse);
	}

out:
	*gamma = exp((log_W0 + (max_eta0 - *scale)) - *log_W);
	assert(*gamma >= 0.0);
	assert(isfinite(*log_W));
}

static void cohort_model_set(struct cohort_model *cm,
			     const struct design *design,
			     const struct vector *recv_coefs)
{
	assert(cm);
	assert(design);
	assert(0 <= cm->isend0 && cm->isend0 < design_send_count(design));
	assert(recv_coefs);
	assert(vector_dim(recv_coefs) == design_recv_dim(design));
	
	ssize_t nreceiver = design_recv_count(design);
	ssize_t dim = design_recv_dim(design);
	ssize_t isend = cm->isend0;
	
	/* The probabilities are p[i] = w[i] / sum(w[j]), so
	 * log(p[i]) = log(w[i]) - log(sum(w[j])).
	 */
	
	/* NOTE: should we worry about overflow in the multiplication? 
	 * It is possible to gaurd against said overflow by computing a
	 * scaled version of eta0: scale the coefficient vector before the
	 * multiplication, then unscale when computing p0.  This
	 * shouldn't be necessary in most (all?) real-world situations.
	 */
	
	/* eta0 */
	design_recv_mul0(1.0, TRANS_NOTRANS, design, isend, recv_coefs, 0.0,
			 &cm->eta0);
	assert(isfinite(vector_norm(&cm->eta0)));
	
	/* max_eta0 */
	cm->max_eta0 = vector_max(&cm->eta0);


	/* store log_p0 in p0 */
	vector_assign_copy(&cm->p0, &cm->eta0);	
	vector_shift(&cm->p0, -cm->max_eta0); /* guard against overflow */
	
	/* log_W0 */
	cm->log_W0 = vector_log_sum_exp(&cm->p0);
	vector_shift(&cm->p0, -cm->log_W0);
	vector_exp(&cm->p0);
	
	/* mean0 */
	design_recv_mul0(1.0, TRANS_TRANS, design, isend, &cm->p0, 0.0,
			 &cm->mean0);
	
	/* imat0 */
	struct vector y;
	struct svector ej;
	double pj;
	ssize_t jrecv;
	
	vector_init(&y, dim);
	svector_init(&ej, nreceiver);
	matrix_fill(&cm->imat0, 0.0);
	
	for (jrecv = 0; jrecv < nreceiver; jrecv++) {
		vector_assign_copy(&y, &cm->mean0);
		svector_set_basis(&ej, jrecv);
		design_recv_muls0(1.0, TRANS_TRANS, design, isend, &ej, -1.0,
				  &y);
		pj = vector_item(&cm->p0, jrecv);
		
		matrix_update1(&cm->imat0, pj, &y, &y);
	}
	svector_deinit(&ej);
	vector_deinit(&y);
	
#ifndef NDEBUG
	/* W0 */
	cm->W0 = exp(cm->log_W0 + cm->max_eta0);
	
	/* w0 */
	vector_assign_copy(&cm->w0, &cm->p0);
	vector_scale(&cm->w0, cm->W0);
#endif
}

static void cohort_model_init(struct cohort_model *cm,
			      const struct design *design,
			      ssize_t isend, const struct vector *recv_coefs)
{
	assert(cm);
	assert(design);
	assert(0 <= isend && isend < design_send_count(design));
	assert(recv_coefs);

	ssize_t nreceiver = design_recv_count(design);
	ssize_t dim = design_recv_dim(design);

	cm->isend0 = isend;
	vector_init(&cm->eta0, nreceiver);		
	vector_init(&cm->p0, nreceiver);
	vector_init(&cm->mean0, dim);
	matrix_init(&cm->imat0, dim, dim);
#ifndef NDEBUG
	vector_init(&cm->w0, nreceiver);
#endif
	cohort_model_set(cm, design, recv_coefs);	
}

static void cohort_model_deinit(struct cohort_model *cm)
{
	assert(cm);
#ifndef NDEBUG
	vector_deinit(&cm->w0);
#endif
	matrix_deinit(&cm->imat0);
	vector_deinit(&cm->mean0);
	vector_deinit(&cm->p0);	
	vector_deinit(&cm->eta0);	
}

static void cohort_models_set(struct intmap *cms,
			      const struct design *design,
			      const struct vector *recv_coefs)
{
	struct intmap_iter it;
	INTMAP_FOREACH(it, cms) {
		struct cohort_model *cm = INTMAP_VAL(it);
		cohort_model_set(cm, design, recv_coefs);
	}
}

static void cohort_models_init(struct intmap *cms,
			       const struct design *design,
			       const struct vector *recv_coefs)
{
	const struct actors *senders = design_senders(design);
	const struct hashset *cohorts = &senders->cohorts;	

	intmap_init(cms, sizeof(struct cohort_model),
		    alignof(struct cohort_model));
	
	struct hashset_iter it;
	HASHSET_FOREACH(it, cohorts) {
		const struct cohort *c = *(struct cohort **)HASHSET_KEY(it);
		struct cohort_iter cit = cohort_iter_make(c);
		if (!cohort_iter_advance(&cit))
			continue; // ignore empty cohorts

		struct cohort_model *cm = intmap_add(cms, (intptr_t)c, NULL);		
		ssize_t isend = COHORT_KEY(cit);
		cohort_model_init(cm, design, isend, recv_coefs);
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
	assert(0 <= isend && isend < model_send_count(model));

	const struct frame *frame = model_frame(model);
	const struct design *design = frame_design(frame);
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
	const struct frame *f = model_frame(m);
	const struct design *d = frame_design(f);

	double max_eta0 = rm->cohort->max_eta0;
	double log_W0 = rm->cohort->log_W0;

	svector_clear(&rm->deta);
	array_clear(&rm->active);

	/* take the initial values if there are self-loops */
	if (design_loops(d)) {
		rm->gamma = 1.0;
		rm->log_W = log_W0;
		rm->scale = max_eta0;
	} else {
		/* add the loop to the active set */
		array_add(&rm->active, &isend);

		/* compute the eta values */
		svector_set_item(&rm->deta, rm->isend, -INFINITY);
		
		/* compute the changes in weights */	
		const struct vector *recv_coefs = model_recv_coefs(rm->model);
		compute_weight_changes(f, isend, recv_coefs, &rm->cohort->eta0,
				       max_eta0, log_W0,
				       &rm->deta, &rm->scale, &rm->gamma, &rm->log_W);

	}
}

static void recv_model_set(struct recv_model *rm,
			   const struct frame *f, 
			   const struct vector *recv_coefs)
{
	recv_model_clear(rm);
	
	const struct design *design = model_design(rm->model);
	bool has_loops = design_loops(design);
	ssize_t dyn_off = design_recv_dyn_index(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	const struct vector beta = vector_slice(recv_coefs, dyn_off, dyn_dim);
	
	/* compute the eta values */
	frame_recv_dmul(1.0, TRANS_NOTRANS, f, rm->isend, &beta, 0.0, &rm->deta);
	if (!has_loops) {
		svector_set_item(&rm->deta, rm->isend, -INFINITY);
	}
	
	/* compute the changes in weights */	
	compute_weight_changes(f, rm->isend, recv_coefs, &rm->cohort->eta0,
			       rm->cohort->max_eta0, rm->cohort->log_W0,
			       &rm->deta, &rm->scale, &rm->gamma, &rm->log_W);

	assert(rm->gamma >= 0.0);
	assert(isfinite(rm->log_W));	

	/* compute the active set */
	struct svector_iter it;
	SVECTOR_FOREACH(it, &rm->deta) {
		ssize_t ix = SVECTOR_IDX(it);
		ssize_t k = array_binary_search(&rm->active, &ix, ssize_compare);
		if (k < 0) {
			array_insert(&rm->active, ~k, &ix);
		}
	}
	assert(array_count(&rm->active) == svector_count(&rm->deta));
}

static void recv_model_init(struct recv_model *rm, struct model *m,
			    ssize_t isend)
{
	assert(rm);
	assert(m);
	assert(0 <= isend && isend < model_send_count(m));

	rm->model = m;
	rm->isend = isend;
	rm->cohort = model_cohort_model(m, isend);
	svector_init(&rm->deta, model_recv_count(m));
	array_init(&rm->active, sizeof(ssize_t));
	recv_model_clear(rm);
}

static void recv_model_deinit(struct recv_model *rm)
{
	assert(rm);
	array_deinit(&rm->active);
	svector_deinit(&rm->deta);
}

static struct recv_model *model_recv_model_raw(struct model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < model_send_count(m));

	struct intmap_pos pos;
	struct recv_model *rm;

	if (!(rm = intmap_find(&m->recv_models, isend, &pos))) {
		rm = intmap_insert(&m->recv_models, &pos, NULL);
		recv_model_init(rm, m, isend);
	}

	assert(rm);
	return rm;
}


static void process_recv_var_event(struct model *m, const struct frame *f, const struct frame_event *e)
{
	assert(m);
	assert(e->type == RECV_VAR_EVENT);
	
	const struct recv_var_event_meta *meta = &e->meta.recv_var;
	struct recv_model *rm = model_recv_model_raw(m, meta->item.isend);
	const struct cohort_model *cm = rm->cohort;
	const struct vector *coefs = model_recv_coefs(m);
	
	//if (rm->isend == 119) {
	//	printf("PROCESS\n");
	//}
	
	ssize_t jrecv = meta->item.jrecv;
	struct svector_pos pos;
	double *pdeta = svector_find(&rm->deta, jrecv, &pos);
	
	if (!pdeta) {
		pdeta = svector_insert(&rm->deta, &pos, 0.0);
		ssize_t ix =
		array_binary_search(&rm->active, &jrecv, ssize_compare);
		assert(ix < 0);
		array_insert(&rm->active, ~ix, &jrecv);
	}
	
	ssize_t index = meta->index;
	double dx = meta->delta;
	double deta = vector_item(coefs, index) * dx;
	
	double scale = rm->scale;
	double gamma = rm->gamma;
	double log_W = rm->log_W;
	double eta = vector_item(&cm->eta0, jrecv) + (*pdeta);
	double eta1 = eta + deta;
	
	/* W' = W + exp(eta') - exp(eta)

	   log(W') = log{ W * (1 + [exp(eta') - exp(eta)]/W) }
		   = log(W) + log[ 1 + (exp(eta') - exp(eta))/W ]
	 
	   gamma' = W0 / W * W / W'
	          = gamma / { 1 + [exp(eta') - exp(eta)] / W }
		  = gamma / { 1 + [ exp(eta' - log(W)) - exp(eta - log(W)) ] }
	 */
	double w = exp(eta - scale);
	double w1 = exp(eta1 - scale);
	double dw = (w1 - w) / exp(log_W);
	double gamma1 = gamma / (1.0 + dw);
	double log_W1 = log_W + log1p(dw);
	
	//if (rm->isend == 119) {
	//	printf("dw: %.22f, eta1: %.22f max_eta: %.22f log_W1: %.22f\n", dw, eta1, max_eta, log_W1);
	//}
	
	/* for dramatic changes in the weight sums, we recompute everything */
	if (!(fabs(dw) <= 0.5)) {
		/* Recompute the diffs when there is overflow */
		//fprintf(stderr, "."); fflush(stderr);
		recv_model_set(rm, f, model_recv_coefs(m));
	} else {
		rm->gamma = gamma1;
		rm->log_W = log_W1;
		assert(rm->gamma >= 0);
		assert(isfinite(log_W1));
		*pdeta += deta;
	}
}


static void handle_frame_event(void *udata, const struct frame_event *e, struct frame *f)
{
	struct model *m = udata;
	assert(m);
	assert(e);
	assert(f == model_frame(m));
	
	switch (e->type) {
		case RECV_VAR_EVENT:
			process_recv_var_event(m, f, e);
			break;
		default:
			break;		/* pass */
	}
}

static void model_clear(struct model *m)
{
	assert(m);
	
	struct intmap_iter it;
	
	INTMAP_FOREACH(it, &m->recv_models) {
		recv_model_clear(INTMAP_VAL(it));
	}
}

static void handle_frame_clear(void *udata, const struct frame *f)
{
	struct model *m = udata;
	assert(m);
	assert(f == model_frame(m));
	model_clear(m);
}

void model_init(struct model *model, struct frame *f,
		const struct vector *recv_coefs)
{
	assert(model);
	assert(f);
	assert(!recv_coefs || design_recv_dim(frame_design(f)) == vector_dim(recv_coefs));
	assert(design_recv_count(frame_design(f)) > 0);
	assert(!design_loops(frame_design(f)) || design_recv_count(frame_design(f)) > 1);

	const struct design *d = frame_design(f);
	
	model->frame = f;
	
	if (recv_coefs) {
		vector_init_copy(&model->recv_coefs, recv_coefs);
	} else {
		vector_init(&model->recv_coefs, design_recv_dim(d));
	}
	
	cohort_models_init(&model->cohort_models, d, &model->recv_coefs);
	intmap_init(&model->recv_models, sizeof(struct recv_model),
		    alignof(struct recv_model));
	refcount_init(&model->refcount);
	
	struct frame_handlers h;
	h.event_mask = RECV_VAR_EVENT;
	h.handle_event = handle_frame_event;
	h.handle_clear = handle_frame_clear;
	frame_add_observer(f, model, &h);
	
}

void model_deinit(struct model *model)
{
	assert(model);

	frame_remove_observer(model->frame, model);
	
	refcount_deinit(&model->refcount);

	struct intmap_iter it;
	INTMAP_FOREACH(it, &model->recv_models) {
		recv_model_deinit(INTMAP_VAL(it));
	}

	intmap_deinit(&model->recv_models);
	cohort_models_deinit(&model->cohort_models);
	vector_deinit(&model->recv_coefs);
}

struct model *model_alloc(struct frame *f, const struct vector *coefs)
{
	assert(f);
	assert(coefs);

	struct model *model = xcalloc(1, sizeof(*model));
	model_init(model, f, coefs);
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

const struct frame *model_frame(const struct model *model)
{
	assert(model);
	return model->frame;
}

const struct design *model_design(const struct model *model)
{
	return frame_design(model_frame(model));
}

const struct vector *model_recv_coefs(const struct model *model)
{
	assert(model);
	return &((struct model *)model)->recv_coefs;
}

ssize_t model_send_count(const struct model *model)
{
	assert(model);
	const struct design *design = model_design(model);
	return design_send_count(design);
}

ssize_t model_recv_count(const struct model *model)
{
	assert(model);
	const struct design *design = model_design(model);
	return design_recv_count(design);
}

ssize_t model_recv_dim(const struct model *model)
{
	assert(model);
	const struct design *design = model_design(model);
	return design_recv_dim(design);
}

struct vector *recv_model_logweight0(const struct recv_model *model)
{
	assert(model);
	return &model->cohort->eta0;
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


void model_set_recv_coefs(struct model *m, const struct vector *recv_coefs)
{
	assert(m);
	assert(!recv_coefs || design_recv_dim(model_design(m)) == vector_dim(recv_coefs));

	if (recv_coefs) {
		vector_assign_copy(&m->recv_coefs, recv_coefs);
	} else {
		vector_fill(&m->recv_coefs, 0.0);
	}
	
	const struct frame *f = model_frame(m);
	const struct design *d = frame_design(f);

	cohort_models_set(&m->cohort_models, d, &m->recv_coefs);
	
	struct intmap_iter it;
	
	INTMAP_FOREACH(it, &m->recv_models) {
		recv_model_set(INTMAP_VAL(it), f, &m->recv_coefs);
	}

}

struct recv_model *model_recv_model(const struct model *m, ssize_t isend)
{
	assert(m);

	struct recv_model *rm = model_recv_model_raw((struct model *)m, isend);
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
 * log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
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

	double scale = rm->scale;
	double log_W = rm->log_W;
	double eta0 = vector_item(&rm->cohort->eta0, jrecv);
	double deta = svector_item(&rm->deta, jrecv);
	double eta = eta0 + deta;
	double log_p = (eta - scale) - log_W;
	
	log_p = MIN(0.0, log_p);
	
	assert(log_p <= 0.0);
	return log_p;
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
			   const struct recv_model *rm, struct vector *y)
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
