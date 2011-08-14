#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "logsumexp.h"
#include "recv_model.h"
#include "util.h"

DEFINE_COMPARE_FN(ssize_compare, ssize_t)

static void
compute_weight_changes(const struct frame *f,
		       const struct recv_model_cohort *cm,
		       const struct svector *deta, double *scale,
		       double *gamma, double *log_W, double *pW)
{
	const struct design *design = f->design;
	const struct vector *eta0 = &cm->eta0;
	double max_eta0 = cm->max_eta0;
	double log_W0 = cm->log_W0;
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

	// fast:
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
			*pW = W + 1;
			*log_W = log1p(W);
		} else {
			*pW = W;
			*log_W = log(W);
		}

		if (*log_W > 0)
			goto out;
	}

	// accurate:
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
		*pW = exp(*log_W);
	}

out:
	*gamma = exp((log_W0 + (max_eta0 - *scale)) - *log_W);
	assert(*gamma >= 0.0);
	assert(isfinite(*log_W));
}

static void cohort_set(struct recv_model_cohort *cm,
		       const struct design *design,
		       const struct vector *recv_coefs)
{
	assert(cm);
	assert(design);
	assert(recv_coefs);
	assert(vector_dim(recv_coefs) == design_recv_dim(design));

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

	/* eta0 */
	design_recv_mul0(1.0, TRANS_NOTRANS, design, recv_coefs, 0.0,
			 &cm->eta0);
	assert(isfinite(vector_max_abs(&cm->eta0)));

	/* max_eta0 */
	cm->max_eta0 = vector_max(&cm->eta0);

	/* store log_p0 in p0 */
	vector_assign_copy(&cm->p0, &cm->eta0);
	vector_shift(&cm->p0, -cm->max_eta0);	/* guard against overflow */

	/* log_W0 */
	cm->log_W0 = vector_log_sum_exp(&cm->p0);
	vector_shift(&cm->p0, -cm->log_W0);
	vector_exp(&cm->p0);

	/* mean0 */
	design_recv_mul0(1.0, TRANS_TRANS, design, &cm->p0, 0.0, &cm->mean0);

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
		design_recv_muls0(1.0, TRANS_TRANS, design, &ej, -1.0, &y);
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

static void cohort_init(struct recv_model_cohort *cm,
			const struct design *design,
			const struct vector *recv_coefs)
{
	assert(cm);
	assert(design);
	assert(recv_coefs);

	ssize_t nreceiver = design_recv_count(design);
	ssize_t dim = design_recv_dim(design);

	vector_init(&cm->eta0, nreceiver);
	vector_init(&cm->p0, nreceiver);
	vector_init(&cm->mean0, dim);
	matrix_init(&cm->imat0, dim, dim);
#ifndef NDEBUG
	vector_init(&cm->w0, nreceiver);
#endif
	cohort_set(cm, design, recv_coefs);
}

static void cohort_deinit(struct recv_model_cohort *cm)
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

ssize_t recv_model_cohort(const struct recv_model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	const struct actor *actors = actors_items(m->senders);
	ssize_t c = actors[isend].cohort;
	return c;
}

static void sender_clear(const struct recv_model *m,
			 struct recv_model_sender *send, ssize_t isend)
{
	assert(m);
	assert(send);
	assert(0 <= isend && isend < recv_model_send_count(m));

	const struct frame *f = recv_model_frame(m);
	const struct design *d = frame_design(f);
	ssize_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	double max_eta0 = cm->max_eta0;
	double log_W0 = cm->log_W0;
	double W0 = cm->W0;

	svector_clear(&send->deta);
	array_clear(&send->active);

	/* take the initial values if there are self-loops */
	if (design_loops(d)) {
		send->gamma = 1.0;
		send->scale = max_eta0;		
		send->log_W = log_W0;
		send->W = W0;
	} else {
		/* add the loop to the active set */
		array_add(&send->active, &isend);

		/* compute the eta values */
		svector_set_item(&send->deta, isend, -INFINITY);

		/* compute the changes in weights */
		compute_weight_changes(f, cm,
				       &send->deta, &send->scale, &send->gamma,
				       &send->log_W, &send->W);

	}
}

static void sender_set(const struct recv_model *m,
		       struct recv_model_sender *send,
		       ssize_t isend,
		       const struct frame *f, const struct vector *recv_coefs)
{
	sender_clear(m, send, isend);

	const struct design *design = recv_model_design(m);
	bool has_loops = design_loops(design);
	ssize_t dyn_off = design_recv_dyn_index(design);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	const struct vector beta = vector_slice(recv_coefs, dyn_off, dyn_dim);

	/* compute the eta values */
	frame_recv_dmul(1.0, TRANS_NOTRANS, f, isend, &beta, 0.0, &send->deta);
	if (!has_loops) {
		svector_set_item(&send->deta, isend, -INFINITY);
	}

	/* compute the changes in weights */
	ssize_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	compute_weight_changes(f, cm,
			       &send->deta, &send->scale, &send->gamma,
			       &send->log_W, &send->W);

	assert(send->gamma >= 0.0);
	assert(isfinite(send->log_W));

	/* compute the active set */
	struct svector_iter it;
	SVECTOR_FOREACH(it, &send->deta) {
		ssize_t ix = SVECTOR_IDX(it);
		ssize_t k =
		    array_binary_search(&send->active, &ix, ssize_compare);
		if (k < 0) {
			array_insert(&send->active, ~k, &ix);
		}
	}
	assert(array_count(&send->active) == svector_count(&send->deta));
}

static void sender_init(const struct recv_model *m,
			struct recv_model_sender *send, ssize_t isend)
{
	assert(send);
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));

	svector_init(&send->deta, recv_model_count(m));
	array_init(&send->active, sizeof(ssize_t));
	sender_clear(m, send, isend);
}

static void sender_deinit(struct recv_model_sender *rm)
{
	assert(rm);
	array_deinit(&rm->active);
	svector_deinit(&rm->deta);
}

static struct recv_model_sender *sender_raw(struct recv_model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));

	return &m->sender_models[isend];
}

static void recv_model_recv_update(void *udata, struct frame *f, ssize_t isend,
				   ssize_t jrecv, ssize_t dyn_index, double delta)
{
	struct recv_model *m = udata;
	const struct actor *actors = actors_items(m->senders);
	ssize_t icohort = actors[isend].cohort;
	struct recv_model_sender *send = sender_raw(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[icohort];
	const struct vector coefs = matrix_col(recv_model_coefs(m), icohort);

	double *pdeta = svector_item_ptr(&send->deta, jrecv);
	//double *pdeta = svector_find(&send->deta, jrecv, &pos);

	if (svector_count(&send->deta) != array_count(&send->active)) {
		assert(svector_count(&send->deta) == array_count(&send->active) + 1);
		ssize_t ix = pdeta - svector_data_ptr(&send->deta);
		assert(ix >= 0);
		assert(svector_index_ptr(&send->deta)[ix] == jrecv);
		array_insert(&send->active, ix, &jrecv);
	}

	const struct design *d = recv_model_design(m);
	ssize_t dyn_off = design_recv_dyn_index(d);
	ssize_t index = dyn_off + dyn_index;

	double deta = vector_item(&coefs, index) * delta;
	double scale = send->scale;
	double gamma = send->gamma;
	double log_W = send->log_W;
	double W = send->W;	
	double eta = vector_item(&cm->eta0, jrecv) + (*pdeta);
	double eta1 = eta + deta;

	/* W' = W + exp(eta') - exp(eta)
	 *
	 * log(W') = log{ W * (1 + [exp(eta') - exp(eta)]/W) }
	 *         = log(W) + log[ 1 + (exp(eta') - exp(eta))/W ]
	 *
	 * gamma' = W0 / W * W / W'
	 *        = gamma / { 1 + [exp(eta') - exp(eta)] / W }
	 *        = gamma / { 1 + [ exp(eta' - log(W)) - exp(eta - log(W)) ] }
	 */
	double w = exp(eta - scale);
	double w1 = exp(eta1 - scale);
	double dw = (w1 - w) / W;
	double gamma1 = gamma / (1.0 + dw);
	double log_W1 = log_W + log1p(dw);
	double W1 = W * dw + W;

	//if (isend == 119) {
	//      printf("dw: %.22f, eta1: %.22f max_eta: %.22f log_W1: %.22f\n", dw, eta1, max_eta, log_W1);
	//}

	/* for dramatic changes in the weight sums, we recompute everything */
	if (!(fabs(dw) <= 0.5)) {
		/* Recompute the diffs when there is overflow */
		//fprintf(stderr, "."); fflush(stderr);
		sender_set(m, send, isend, f, &coefs);
	} else {
		send->gamma = gamma1;
		send->log_W = log_W1;
		send->W = W1;
		assert(send->gamma >= 0);
		assert(isfinite(log_W1));
		*pdeta += deta;
	}
}

static void model_clear(struct recv_model *m)
{
	assert(m);

	struct recv_model_sender *senders = m->sender_models;
	ssize_t i, n = recv_model_send_count(m);

	for (i = 0; i < n; i++) {
		sender_clear(m, &senders[i], i);
	}
}

static void recv_model_clear(void *udata, struct frame *f)
{
	struct recv_model *m = udata;
	assert(m);
	assert(f == recv_model_frame(m));
	model_clear(m);
}

void recv_model_init(struct recv_model *model, struct frame *f,
		     const struct actors *senders, const struct matrix *coefs)
{
	assert(model);
	assert(f);
	assert(senders);
	assert(actors_count(senders) == design_send_count(frame_design(f)));
	assert(!coefs
	       || design_recv_dim(frame_design(f)) == matrix_nrow(coefs));
	assert(!coefs || actors_cohort_count(senders) == matrix_ncol(coefs));
	assert(design_recv_count(frame_design(f)) > 0);
	assert(!design_loops(frame_design(f))
	       || design_recv_count(frame_design(f)) > 1);

	const struct design *d = frame_design(f);
	ssize_t ncohort = actors_cohort_count(senders);

	model->frame = f;
	model->senders = senders;

	if (coefs) {
		matrix_init_copy(&model->coefs, TRANS_NOTRANS, coefs);
	} else {
		matrix_init(&model->coefs, design_recv_dim(d), ncohort);
	}

	struct recv_model_cohort *cms = xcalloc(ncohort, sizeof(*cms));
	ssize_t ic;
	for (ic = 0; ic < ncohort; ic++) {
		struct vector col = matrix_col(&model->coefs, ic);
		cohort_init(&cms[ic], d, &col);
	}
	model->cohort_models = cms;

	ssize_t isend, nsend = design_send_count(d);
	struct recv_model_sender *sms = xcalloc(nsend, sizeof(*sms));

	for (isend = 0; isend < nsend; isend++) {
		sender_init(model, &sms[isend], isend);
	}
	model->sender_models = sms;

	struct frame_callbacks callbacks = {
		NULL,		// msg_add
		NULL,		// msg_advance
		recv_model_recv_update,
		NULL,		// send_update
		recv_model_clear
	};

	frame_add_observer(f, model, &callbacks);
}

void recv_model_deinit(struct recv_model *model)
{
	assert(model);

	frame_remove_observer(model->frame, model);

	struct recv_model_sender *sms = model->sender_models;
	ssize_t isend, nsend = recv_model_send_count(model);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&sms[isend]);
	}
	xfree(sms);

	struct recv_model_cohort *cms = model->cohort_models;
	ssize_t ic, nc = recv_model_cohort_count(model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cms[ic]);
	}
	xfree(cms);

	matrix_deinit(&model->coefs);
}

const struct frame *recv_model_frame(const struct recv_model *model)
{
	assert(model);
	return model->frame;
}

const struct design *recv_model_design(const struct recv_model *model)
{
	return frame_design(recv_model_frame(model));
}

const struct matrix *recv_model_coefs(const struct recv_model *model)
{
	assert(model);
	return &((struct recv_model *)model)->coefs;
}

ssize_t recv_model_send_count(const struct recv_model *model)
{
	assert(model);
	const struct design *design = recv_model_design(model);
	return design_send_count(design);
}

ssize_t recv_model_cohort_count(const struct recv_model *model)
{
	assert(model);
	return actors_cohort_count(model->senders);
}

ssize_t recv_model_count(const struct recv_model *model)
{
	assert(model);
	const struct design *design = recv_model_design(model);
	return design_recv_count(design);
}

ssize_t recv_model_dim(const struct recv_model *model)
{
	assert(model);
	const struct design *design = recv_model_design(model);
	return design_recv_dim(design);
}

double recv_model_logsumwt0(const struct recv_model *m, ssize_t c)
{
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));
	return m->cohort_models[c].log_W0 + m->cohort_models[c].max_eta0;
}

struct vector *recv_model_logwts0(const struct recv_model *m, ssize_t c)
{
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));
	return &((struct recv_model *)m)->cohort_models[c].eta0;
}

struct vector *recv_model_probs0(const struct recv_model *m, ssize_t c)
{
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));
	return &((struct recv_model *)m)->cohort_models[c].p0;
}

double recv_model_prob0(const struct recv_model *m, ssize_t c, ssize_t jrecv)
{
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));
	assert(0 <= jrecv && jrecv < recv_model_count(m));

	const struct vector *p0 = &m->cohort_models[c].p0;
	return vector_item(p0, jrecv);
}

struct vector *recv_model_mean0(const struct recv_model *m, ssize_t c)
{
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));
	return &((struct recv_model *)m)->cohort_models[c].mean0;
}

struct matrix *recv_model_imat0(const struct recv_model *m, ssize_t c)
{
	assert(m);
	assert(0 <= c && c < recv_model_cohort_count(m));
	return &((struct recv_model *)m)->cohort_models[c].imat0;
}

void recv_model_set_coefs(struct recv_model *m, const struct matrix *coefs)
{
	assert(m);
	assert(!coefs
	       || design_recv_dim(recv_model_design(m)) == matrix_nrow(coefs));

	if (coefs) {
		matrix_assign_copy(&m->coefs, TRANS_NOTRANS, coefs);
	} else {
		matrix_fill(&m->coefs, 0.0);
	}

	const struct frame *f = recv_model_frame(m);
	const struct design *d = frame_design(f);
	const struct cohort *cohorts = actors_cohorts(m->senders);
	struct recv_model_sender *senders = m->sender_models;

	ssize_t ic, nc = recv_model_cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		struct vector col = matrix_col(&m->coefs, ic);
		cohort_set(&m->cohort_models[ic], d, &col);

		const ssize_t *isend = array_to_ptr(&cohorts[ic].actors);
		ssize_t i, n = array_count(&cohorts[ic].actors);
		for (i = 0; i < n; i++) {
			sender_set(m, &senders[isend[i]], isend[i], f, &col);
		}
	}
}

struct recv_model_sender *recv_model_send(const struct recv_model *m,
					  ssize_t isend)
{
	assert(m);

	struct recv_model_sender *rm =
	    sender_raw((struct recv_model *)m, isend);
	return rm;
}

void recv_model_get_active(const struct recv_model *m, ssize_t isend,
			   ssize_t **jrecv, ssize_t *n)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	assert(jrecv);
	assert(n);

	const struct recv_model_sender *rm = recv_model_send(m, isend);
	*jrecv = array_to_ptr(&rm->active);
	*n = array_count(&rm->active);
}

double recv_model_logsumwt(const struct recv_model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));

	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->log_W + rm->scale;
}

double recv_model_invgrow(const struct recv_model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));

	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->gamma;
}

/*
 * log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double recv_model_logprob(const struct recv_model *m, ssize_t isend,
			  ssize_t jrecv)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	assert(0 <= jrecv && jrecv < recv_model_count(m));

	ssize_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	const struct recv_model_sender *send = recv_model_send(m, isend);
	/*
	   double gamma = ctx->gamma;
	   double p0 = vector_item(ctx->group->p0, jrecv);
	   double dp = svector_item(ctx->dp, jrecv);
	   double p = gamma * p0 + dp;
	   p = MAX(0.0, MIN(1.0, p));
	   return log(p);
	 */

	double scale = send->scale;
	double log_W = send->log_W;
	double eta0 = vector_item(&cm->eta0, jrecv);
	double deta = svector_item(&send->deta, jrecv);
	double eta = eta0 + deta;
	double log_p = (eta - scale) - log_W;

	log_p = MIN(0.0, log_p);

	assert(log_p <= 0.0);
	return log_p;
}

double recv_model_prob(const struct recv_model *m, ssize_t isend, ssize_t jrecv)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	assert(0 <= jrecv && jrecv < recv_model_count(m));

	double lp = recv_model_logprob(m, isend, jrecv);
	double p = exp(lp);
	return p;
}

void recv_model_axpy_probs(double alpha,
			   const struct recv_model *m,
			   ssize_t isend, struct vector *y)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	assert(y);
	assert(vector_dim(y) == recv_model_count(m));

	ssize_t j, n = recv_model_count(m);

	for (j = 0; j < n; j++) {
		double p = recv_model_prob(m, isend, j);
		*vector_item_ptr(y, j) += alpha * p;
	}
}

const struct actors *recv_model_senders(const struct recv_model *model)
{
	assert(model);
	return model->senders;
}
