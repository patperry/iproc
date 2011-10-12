#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "logsumexp.h"
#include "sblas.h"
#include "recv_model.h"


static void
compute_weight_changes(struct recv_model_sender *sm,
		       const struct recv_model_cohort *cm)
{
	const struct vector *eta0 = &cm->eta0;
	double max_eta0 = cm->max_eta0;
	double log_W0 = cm->log_W0;
	ssize_t jrecv, nrecv = vector_dim(eta0);
	bool shrink = false;
	size_t ia, na = sm->nactive;

	/* compute the maximum eta value */
	sm->scale = max_eta0;

	for (ia = 0; ia < na; ia++) {
		ptrdiff_t jrecv = sm->active[ia];
		double eta0_j = vector_item(eta0, jrecv);
		double deta_j = sm->deta[ia];
		double eta_j = eta0_j + deta_j;

		if (eta0_j == max_eta0 && deta_j < 0)
			shrink = true;

		sm->scale = MAX(sm->scale, eta_j);
	}

	// fast:
	{
		double W = exp(log_W0 + (max_eta0 - sm->scale));
		bool found_max = false;

		for (ia = 0; ia < na; ia++) {
			ssize_t jrecv = sm->active[ia];
			double eta0_j = vector_item(eta0, jrecv);
			double deta_j = sm->deta[ia];
			double eta_j = eta0_j + deta_j;

			if (!found_max && eta_j == sm->scale) {
				found_max = true;
				W += -exp(eta0_j - sm->scale);
			} else {
				W += (exp(eta_j - sm->scale)
				      - exp(eta0_j - sm->scale));
			}
		}

		if (found_max) {
			sm->W = W + 1;
			sm->log_W = log1p(W);
		} else {
			sm->W = W;
			sm->log_W = log(W);
		}

		if (sm->log_W > 0)
			goto out;
	}

	// accurate:
	{
		//fprintf(stderr, "!"); fflush(stderr);
		/* compute the new max_eta */
		ia = 0;
		if (shrink && sm->scale == max_eta0) {
			double max = -INFINITY;
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double eta0_j = vector_item(eta0, jrecv);
				double deta_j = 0.0;

				if (ia < na && sm->active[ia] == jrecv) {
					deta_j = sm->deta[ia];
					ia++;
				}		
				
				double eta_j = eta0_j + deta_j;

				max = MAX(max, eta_j);
			}
			sm->scale = max;
		}

		/* compute the new log_W */
		struct logsumexp lse;
		logsumexp_init(&lse);

		ia = 0;
		for (jrecv = 0; jrecv < nrecv; jrecv++) {
			double eta0_j = vector_item(eta0, jrecv);
			double deta_j = 0.0;
		
			if (ia < na && sm->active[ia] == jrecv) {
				deta_j = sm->deta[ia];
				ia++;
			}	

			double eta_j = eta0_j + deta_j;

			logsumexp_insert(&lse, eta_j - sm->scale);
		}
		sm->log_W = logsumexp_value(&lse);
		sm->W = exp(sm->log_W);
	}

out:
	sm->gamma = exp((log_W0 + (max_eta0 - sm->scale)) - sm->log_W);
	assert(sm->gamma >= 0.0);
	assert(isfinite(sm->log_W));
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
	design_recv_mul0(1.0, BLAS_NOTRANS, design, recv_coefs, 0.0,
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
	design_recv_mul0(1.0, BLAS_TRANS, design, &cm->p0, 0.0, &cm->mean0);

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
		design_recv_muls0(1.0, BLAS_TRANS, design, &ej, -1.0, &y);
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

ptrdiff_t recv_model_cohort(const struct recv_model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	return m->cohorts[isend];
}

static void sender_grow_active(struct recv_model_sender *send)
{
	size_t nactive_max = send->nactive_max;

	if (send->nactive == nactive_max) {
		nactive_max = ARRAY_GROW(nactive_max, SIZE_MAX);
		send->deta = xrealloc(send->deta, nactive_max *
				sizeof(send->deta[0]));
		send->active = xrealloc(send->active, nactive_max *
				sizeof(send->active[0]));
		send->nactive_max = nactive_max;
	}
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

	send->nactive = 0;

	/* take the initial values if there are self-loops */
	if (design_loops(d)) {
		send->gamma = 1.0;
		send->scale = max_eta0;		
		send->log_W = log_W0;
		send->W = W0;
	} else {
		sender_grow_active(send);

		/* add the loop to the active set */
		send->active[0] = isend;
		send->deta[0] = -INFINITY;
		send->nactive = 1;

		/* compute the changes in weights */
		compute_weight_changes(send, cm);
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
	struct vector *dx;
	ptrdiff_t *active;
	size_t iz, nz;
	frame_recv_get_dx(f, isend, &dx, &active, &nz);

	send->nactive = 0;
	for (iz = 0; iz < nz; iz++) {
		sender_grow_active(send);
		send->active[iz] = active[iz];
		send->deta[iz] = vector_dot(&beta, &dx[iz]);
		send->nactive++;
	}

	if (!has_loops) {
		ptrdiff_t ix = sblas_find(send->nactive, send->active, isend);
		if (ix < 0) {
			sender_grow_active(send);
			ix = ~ix;
			memmove(send->active + ix + 1, send->active + ix,
				(send->nactive - ix) * sizeof(send->active[0]));
			memmove(send->deta + ix + 1, send->deta + ix,
				(send->nactive - ix) * sizeof(send->deta[0]));
			send->active[ix] = isend;
			send->nactive++;
		}
		send->deta[ix] = -INFINITY;
	}

	/* compute the changes in weights */
	ssize_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	compute_weight_changes(send, cm);

	assert(send->gamma >= 0.0);
	assert(isfinite(send->log_W));
}

static void sender_init(const struct recv_model *m,
			struct recv_model_sender *send, ssize_t isend)
{
	assert(send);
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));

	send->deta = NULL;
	send->active = NULL;
	send->nactive = 0;
	send->nactive_max = 0;
	sender_clear(m, send, isend);
}

static void sender_deinit(struct recv_model_sender *send)
{
	assert(send);
	free(send->active);
	free(send->deta);
}

static struct recv_model_sender *sender_raw(struct recv_model *m, ssize_t isend)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));

	return &m->sender_models[isend];
}

static void recv_model_recv_update(void *udata, struct frame *f, ssize_t isend,
				   ssize_t jrecv, const struct svector *delta)
{
	struct recv_model *m = udata;
	const struct design *d = frame_design(f);
	assert(recv_model_design(m) == d);
	
	if (isend == jrecv && !design_loops(d))
		return;

	ptrdiff_t icohort = recv_model_cohort(m, isend);
	struct recv_model_sender *send = sender_raw(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[icohort];
	const struct vector coefs = matrix_col(recv_model_coefs(m), icohort);

	ptrdiff_t ix = sblas_find(send->nactive, send->active, jrecv);
	if (ix < 0) {
		sender_grow_active(send);
		ix = ~ix;
		memmove(send->active + ix + 1, send->active + ix,
			(send->nactive - ix) * sizeof(send->active[0]));
		memmove(send->deta + ix + 1, send->deta + ix,
			(send->nactive - ix) * sizeof(send->deta[0]));
		send->active[ix] = jrecv;
		send->deta[ix] = 0.0;
		send->nactive++;
	}
	double *pdeta = &send->deta[ix];

	ssize_t dyn_off = design_recv_dyn_index(d);
	ssize_t dyn_dim = design_recv_dyn_dim(d);
	struct vector dyn_coefs = vector_slice(&coefs, dyn_off, dyn_dim);

	double deta = svector_dot(delta, &dyn_coefs);
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
	(void)f; // unused
	struct recv_model *m = udata;
	assert(m);
	assert(f == recv_model_frame(m));
	model_clear(m);
}


void recv_model_init(struct recv_model *model,
		     struct frame *f,
		     size_t ncohort,
		     const ptrdiff_t *cohorts,
		     const struct matrix *coefs)
{
	assert(model);
	assert(f);
	assert(ncohort > 0);
	assert(cohorts);
	assert(!coefs
	       || design_recv_dim(frame_design(f)) == matrix_nrow(coefs));
	assert(!coefs || ncohort);
	assert(design_recv_count(frame_design(f)) > 0);
	assert(!design_loops(frame_design(f))
	       || design_recv_count(frame_design(f)) > 1);

	const struct design *d = frame_design(f);
	const size_t nsend = design_send_count(d);

	model->frame = f;
	model->ncohort = ncohort;
	model->cohorts = xmalloc(nsend * sizeof(cohorts[0]));
	memcpy(model->cohorts, cohorts, sizeof(cohorts[0]) * nsend);

	if (coefs) {
		matrix_init_copy(&model->coefs, BLAS_NOTRANS, coefs);
	} else {
		matrix_init(&model->coefs, design_recv_dim(d), ncohort);
	}

	struct recv_model_cohort *cms = xcalloc(ncohort, sizeof(*cms));
	size_t ic;
	for (ic = 0; ic < ncohort; ic++) {
		struct vector col = matrix_col(&model->coefs, ic);
		cohort_init(&cms[ic], d, &col);
	}
	model->cohort_models = cms;

	size_t isend;
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
	free(sms);

	struct recv_model_cohort *cms = model->cohort_models;
	ssize_t ic, nc = recv_model_cohort_count(model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cms[ic]);
	}
	free(cms);

	matrix_deinit(&model->coefs);
	free(model->cohorts);
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
	return model->ncohort;
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
		matrix_assign_copy(&m->coefs, BLAS_NOTRANS, coefs);
	} else {
		matrix_fill(&m->coefs, 0.0);
	}

	const struct frame *f = recv_model_frame(m);
	const struct design *d = frame_design(f);
	struct recv_model_sender *senders = m->sender_models;

	size_t ic, nc = recv_model_cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		struct vector col = matrix_col(&m->coefs, ic);
		cohort_set(&m->cohort_models[ic], d, &col);
	}

	const ptrdiff_t *cohorts = recv_model_cohorts(m);
	size_t isend, nsend = recv_model_send_count(m);
	for (isend = 0; isend < nsend; isend++) {
		ptrdiff_t ic = cohorts[isend];
		struct vector col = matrix_col(&m->coefs, ic);
		sender_set(m, &senders[isend], isend, f, &col);
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
			   ptrdiff_t **jrecv, ssize_t *n)
{
	assert(m);
	assert(0 <= isend && isend < recv_model_send_count(m));
	assert(jrecv);
	assert(n);

	const struct recv_model_sender *sm = recv_model_send(m, isend);
	*jrecv = sm->active;
	*n = sm->nactive;
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

	double deta = 0.0;
	ptrdiff_t iactive = sblas_find(send->nactive, send->active, jrecv);

	if (iactive >= 0) {
		deta = send->deta[iactive];
	}

	double scale = send->scale;
	double log_W = send->log_W;
	double eta0 = vector_item(&cm->eta0, jrecv);
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

const ptrdiff_t *recv_model_cohorts(const struct recv_model *model)
{
	assert(model);
	return model->cohorts;
}

