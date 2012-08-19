#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "logsumexp.h"
#include "lapack.h"
#include "matrixutil.h"
#include "util.h"
#include "sblas.h"
#include "recv_model.h"


void recv_coefs_init(struct recv_coefs *c, const struct frame *f)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_tvar_dim(r);
	size_t dimr = dimr0 + dimr1;	
	size_t dimd0 = design2_trait_dim(d);
	size_t dimd1 = design2_tvar_dim(d);
	size_t dimd = dimd0 + dimd1;
	size_t dim = dimr + dimd;
	double *all = xmalloc(dim * sizeof(*all));
	double *allr = all;
	double *alld = allr + dimr;
	
	c->all = all;
	c->dim = dim;
	
	c->recv.all = allr;
	c->recv.dim = dimr;
	c->recv.traits = allr;
	c->recv.tvars = allr + dimr0;
	
	c->dyad.all = alld;
	c->dyad.dim = dimd;
	c->dyad.traits = alld;
	c->dyad.tvars = alld + dimd0;
}


void recv_coefs_deinit(struct recv_coefs *c)
{
	free(c->all);
}


static void
compute_weight_changes(struct recv_model_sender *sm,
		       const struct recv_model_cohort *cm,
		       size_t nrecv)
{
	const double *eta0 = cm->eta0;
	double max_eta0 = cm->max_eta0;
	double log_W0 = cm->log_W0;
	size_t jrecv;
	int shrink = 0;
	size_t ia, na = sm->active.nz;

	/* compute the maximum eta value */
	sm->scale = max_eta0;

	for (ia = 0; ia < na; ia++) {
		size_t jrecv = sm->active.indx[ia];
		double eta0_j = eta0[jrecv];
		double deta_j = sm->deta[ia];
		double eta_j = eta0_j + deta_j;

		if (eta0_j == max_eta0 && deta_j < 0)
			shrink = 1;

		sm->scale = MAX(sm->scale, eta_j);
	}

	// fast:
	{
		double W = exp(log_W0 + (max_eta0 - sm->scale));
		int found_max = 0;

		for (ia = 0; ia < na; ia++) {
			size_t jrecv = sm->active.indx[ia];
			double eta0_j = eta0[jrecv];
			double deta_j = sm->deta[ia];
			double eta_j = eta0_j + deta_j;

			if (!found_max && eta_j == sm->scale) {
				found_max = 1;
				W += -exp(eta0_j - sm->scale);
			} else {
				W += (exp(eta_j - sm->scale)
				      - exp(eta0_j - sm->scale));
			}
		}

		if (found_max) {
			//sm->W = W + 1;
			sm->log_W = log1p(W);
		} else {
			//sm->W = W;
			sm->log_W = log(W);
		}

		if (sm->log_W > 0)
			goto out;
	}

	// accurate:
	{
		fprintf(stderr, "!"); fflush(stderr);
		/* compute the new max_eta */
		ia = 0;
		if (shrink && sm->scale == max_eta0) {
			double max = -INFINITY;
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double eta0_j = eta0[jrecv];
				double deta_j = 0.0;

				if (ia < na && sm->active.indx[ia] == jrecv) {
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
			double eta0_j = eta0[jrecv];
			double deta_j = 0.0;

			if (ia < na && sm->active.indx[ia] == jrecv) {
				deta_j = sm->deta[ia];
				ia++;
			}

			double eta_j = eta0_j + deta_j;

			logsumexp_insert(&lse, eta_j - sm->scale);
		}
		sm->log_W = logsumexp_value(&lse);
		//sm->W = exp(sm->log_W);
	}

out:
	sm->gamma = exp((log_W0 + (max_eta0 - sm->scale)) - sm->log_W);
	assert(sm->gamma >= 0.0);
	assert(isfinite(sm->log_W));
}

static void cohort_set(struct recv_model_cohort *cm, size_t c,
		       const struct frame *f,
		       const struct recv_coefs *coefs)
{
	const struct design *s = frame_send_design(f);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t isend = design_cohort_rep(s, c);
	size_t nreceiver = frame_recv_count(f);

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
	design_traits_mul(1.0, r, coefs->recv.traits, 0.0, cm->eta0);
	design2_traits_mul(1.0, d, isend, coefs->dyad.traits, 1.0, cm->eta0);	
	
	/* max_eta0 */
	cm->max_eta0 = vector_max(nreceiver, cm->eta0);

	/* store log_p0 in p0 */
	blas_dcopy(nreceiver, cm->eta0, 1, cm->p0, 1);
	vector_shift(nreceiver, -cm->max_eta0, cm->p0);	/* guard against overflow */

	/* log_W0 */
	cm->log_W0 = vector_logsumexp(nreceiver, cm->p0);
	vector_shift(nreceiver, -cm->log_W0, cm->p0);
	vector_exp(nreceiver, cm->p0);

/*
	// mean0
	design_mul0(1.0, BLAS_TRANS, design, cm->p0, 0.0, cm->mean0);

	// imat0
	double one = 1.0;
	double pj;
	size_t jrecv;
	struct vpattern pat_j;
	pat_j.nz = 1;
	pat_j.indx = &jrecv;

	double *y = xmalloc(dim * sizeof(double));
	matrix_dzero(dim, dim, &cm->imat0);

	for (jrecv = 0; jrecv < nreceiver; jrecv++) {
		blas_dcopy(dim, cm->mean0, 1, y, 1);
		design_muls0(1.0, BLAS_TRANS, design, &one, &pat_j, -1.0, y);
		pj = cm->p0[jrecv];

		blas_dger(dim, dim, pj, y, 1, y, 1, &cm->imat0);
	}
	free(y);
*/
 
#ifndef NDEBUG
	/* W0 */
	//cm->W0 = exp(cm->log_W0 + cm->max_eta0);

	/* w0 */
	blas_dcopy(nreceiver, cm->p0, 1, cm->w0, 1);
	blas_dscal(nreceiver, exp(cm->log_W0), cm->w0, 1);
#endif
}

static void cohort_init(struct recv_model_cohort *cm,
			size_t c,
			const struct frame *f,
			const struct recv_coefs *coefs)
{
	size_t nrecv = frame_recv_count(f);
	size_t dim = coefs->recv.dim + coefs->dyad.dim;
	size_t dim2 = dim * dim;

	cm->eta0 = xmalloc(nrecv * sizeof(double));
	cm->p0 = xmalloc(nrecv * sizeof(double));
	cm->mean0 = xmalloc(dim * sizeof(double));
	cm->imat0 = xmalloc(dim2 * sizeof(double));
	cm->ldimat0 = MAX(1, dim);
#ifndef NDEBUG
	cm->w0 = xmalloc(nrecv * sizeof(double));
#endif
	cohort_set(cm, c, f, coefs);
}

static void cohort_deinit(struct recv_model_cohort *cm)
{
	assert(cm);
#ifndef NDEBUG
	free(cm->w0);
#endif
	free(cm->imat0);
	free(cm->mean0);
	free(cm->p0);
	free(cm->eta0);
}

size_t recv_model_cohort(const struct recv_model *m, size_t isend)
{
	const struct frame *f = recv_model_frame(m);
	const struct design *s = frame_send_design(f);
	return design_cohort(s, isend);
}

static void sender_clear(const struct recv_model *m,
			 struct recv_model_sender *send, size_t isend)
{
	assert(m);
	assert(send);
	assert(isend < recv_model_send_count(m));

	size_t nrecv = recv_model_count(m);
	const struct frame *f = recv_model_frame(m);
	size_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	double max_eta0 = cm->max_eta0;
	double log_W0 = cm->log_W0;
	//double W0 = cm->W0;

	vpattern_clear(&send->active);

	/* take the initial values if there are self-loops */
	if (frame_has_loops(f)) {
		send->gamma = 1.0;
		send->scale = max_eta0;
		send->log_W = log_W0;
		//send->W = W0;
	} else {
		size_t nzmax1;
		if ((nzmax1 = vpattern_grow(&send->active, 1))) {
			send->deta = xrealloc(send->deta,
					      nzmax1 * sizeof(send->deta[0]));
		}

		/* add the loop to the active set */
		send->active.indx[0] = isend;
		send->active.nz = 1;
		send->deta[0] = -INFINITY;

		/* compute the changes in weights */
		compute_weight_changes(send, cm, nrecv);
	}
}

static void sender_set(const struct recv_model *m,
		       struct recv_model_sender *send,
		       size_t isend,
		       const struct frame *f, const struct recv_coefs *coefs)
{
	sender_clear(m, send, isend);

	const struct design2 *d = frame_dyad_design(f);
	size_t dyn_dim = design2_tvar_dim(d);
	int has_loops = frame_has_loops(f);
	const double *beta = coefs->dyad.tvars;

	/* compute the eta values */
	const double *dx;
	const size_t *active;
	size_t iz, nz;
	design2_tvars_get(d, isend, &dx, &active, &nz);

	vpattern_clear(&send->active);
	size_t nzmax1;
	if ((nzmax1 = vpattern_grow(&send->active, nz))) {
		send->deta = xrealloc(send->deta,
				      nzmax1 * sizeof(send->deta[0]));
	}

	memcpy(send->active.indx, active, nz * sizeof(active[0]));
	send->active.nz = nz;
	for (iz = 0; iz < nz; iz++) {
		send->deta[iz] = blas_ddot(dyn_dim, beta, 1, dx + iz * dyn_dim, 1);
	}

	if (!has_loops) {
		size_t nzmax = send->active.nzmax;
		int ins;
		size_t ix = vpattern_search(&send->active, isend, &ins);

		if (ins) {
			if (nzmax != send->active.nzmax) {
				nzmax = send->active.nzmax;
				send->deta = xrealloc(send->deta, nzmax * sizeof(send->deta[0]));
			}
			memmove(send->deta + ix + 1, send->deta + ix,
				(send->active.nz - 1 - ix) * sizeof(send->deta[0]));
		}
		send->deta[ix] = -INFINITY;
	}

	/* compute the changes in weights */
	size_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	size_t nrecv = recv_model_count(m);
	compute_weight_changes(send, cm, nrecv);

	assert(send->gamma >= 0.0);
	assert(isfinite(send->log_W));
}

static void sender_init(const struct recv_model *m,
			struct recv_model_sender *send, size_t isend)
{
	assert(send);
	assert(m);
	assert(isend < recv_model_send_count(m));

	send->deta = NULL;
	vpattern_init(&send->active);
	sender_clear(m, send, isend);
}

static void sender_deinit(struct recv_model_sender *send)
{
	assert(send);
	vpattern_deinit(&send->active);
	free(send->deta);
}

static struct recv_model_sender *sender_raw(struct recv_model *m, size_t isend)
{
	assert(m);
	assert(isend < recv_model_send_count(m));

	return &m->sender_models[isend];
}

static void recv_model_dyad_update(void *udata, struct design2 *d, size_t isend,
				   size_t jrecv, const double *delta,
				   const struct vpattern *pat)
{
	struct recv_model *m = udata;
	struct frame *f = design2_frame(d);	

	if (isend == jrecv && !frame_has_loops(f))
		return;

	size_t icohort = recv_model_cohort(m, isend);
	struct recv_model_sender *send = sender_raw(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[icohort];

	size_t nzmax = send->active.nzmax;
	int ins;

	size_t ix = vpattern_search(&send->active, jrecv, &ins);
	if (ins) {
		if (nzmax != send->active.nzmax) {
			nzmax = send->active.nzmax;
			send->deta = xrealloc(send->deta, nzmax *
					      sizeof(send->deta[0]));
		}

		memmove(send->deta + ix + 1, send->deta + ix,
			(send->active.nz - 1 - ix) * sizeof(send->deta[0]));
		send->deta[ix] = 0.0;
	}
	double *pdeta = &send->deta[ix];

	const double *dyn_coefs = m->coefs.dyad.tvars;

	double deta = sblas_ddoti(delta, pat, dyn_coefs);
	double scale = send->scale;
	double gamma = send->gamma;
	double log_W = send->log_W;
	//double W = send->W;
	double eta = cm->eta0[jrecv] + (*pdeta);

	/* W' = W + exp(eta') - exp(eta)
	 *
	 * log(W') = log{ W * (1 + [exp(eta') - exp(eta)]/W) }
	 *         = log(W) + log[ 1 + (exp(eta') - exp(eta))/W ]
	 *
	 * gamma' = W0 / W * W / W'
	 *        = gamma / { 1 + [exp(eta') - exp(eta)] / W }
	 *        = gamma / { 1 + [ exp(eta' - log(W)) - exp(eta - log(W)) ] }
	 */
	double dw = exp(eta - scale - log_W) * expm1(deta);
	double gamma1 = gamma / (1.0 + dw);
	double log_W1 = log_W + log1p(dw);
	//double W1 = W * (1 + dw);	


	//if (isend == 119) {
	//      printf("dw: %.22f, eta1: %.22f max_eta: %.22f log_W1: %.22f\n", dw, eta1, max_eta, log_W1);
	//}

	/* for dramatic changes in the weight sums, we recompute everything */
	if (!(fabs(dw) <= 0.3 * fabs(log_W1))) {
		/* Recompute the diffs when there is overflow */
		fprintf(stderr, "."); fflush(stderr);
		sender_set(m, send, isend, f, &m->coefs);
	} else {
		//fprintf(stderr, "."); fflush(stderr);
		send->gamma = gamma1;
		send->log_W = log_W1;
		//send->W = W1;
		assert(send->gamma >= 0);
		assert(isfinite(log_W1));
		*pdeta += deta;
	}
}

static void model_clear(struct recv_model *m)
{
	assert(m);

	struct recv_model_sender *senders = m->sender_models;
	size_t i, n = recv_model_send_count(m);

	for (i = 0; i < n; i++) {
		sender_clear(m, &senders[i], i);
	}
}

static void recv_model_dyad_clear(void *udata, struct design2 *d)
{
	(void)d;		// unused
	struct recv_model *m = udata;
	assert(m);
	model_clear(m);
}

void recv_model_init(struct recv_model *model,
		     struct frame *f,
		     const struct recv_coefs *coefs)
{
	assert(model);
	assert(f);
	assert(frame_recv_count(f) > 0);
	assert(!frame_has_loops(f) || frame_recv_count(f) > 1);

	const struct design *s = frame_send_design(f);
	struct design2 *d = frame_dyad_design(f);
	const size_t nsend = frame_send_count(f);

	model->frame = f;

	recv_coefs_init(&model->coefs, f);
	
	if (coefs) {
		memcpy(model->coefs.all, coefs->all, model->coefs.dim * sizeof(*model->coefs.all));
	} else {
		memset(model->coefs.all, 0, model->coefs.dim * sizeof(*model->coefs.all));		
	}

	size_t nc = design_cohort_count(s);
	struct recv_model_cohort *cms = xcalloc(nc, sizeof(*cms));
	size_t ic;
	
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cms[ic], ic, f, &model->coefs);
	}
	model->cohort_models = cms;

	size_t isend;
	struct recv_model_sender *sms = xcalloc(nsend, sizeof(*sms));

	for (isend = 0; isend < nsend; isend++) {
		sender_init(model, &sms[isend], isend);
	}
	model->sender_models = sms;

	struct design2_callbacks callbacks = {
		recv_model_dyad_update,
		NULL,
		recv_model_dyad_clear
	};

	design2_add_observer(d, model, &callbacks);
}

void recv_model_deinit(struct recv_model *model)
{
	assert(model);

	frame_remove_observer(model->frame, model);

	struct recv_model_sender *sms = model->sender_models;
	size_t isend, nsend = recv_model_send_count(model);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&sms[isend]);
	}
	free(sms);

	struct recv_model_cohort *cms = model->cohort_models;
	size_t ic, nc = recv_model_cohort_count(model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cms[ic]);
	}
	free(cms);

	recv_coefs_deinit(&model->coefs);
}

struct frame *recv_model_frame(const struct recv_model *model)
{
	assert(model);
	return model->frame;
}

const struct design *recv_model_design(const struct recv_model *model)
{
	return frame_recv_design(recv_model_frame(model));
}

const struct recv_coefs *recv_model_coefs(const struct recv_model *model)
{
	assert(model);
	return &((struct recv_model *)model)->coefs;
}

size_t recv_model_send_count(const struct recv_model *model)
{
	assert(model);
	return frame_send_count(model->frame);
}

size_t recv_model_cohort_count(const struct recv_model *model)
{
	assert(model);
	const struct frame *f = recv_model_frame(model);
	const struct design *s = frame_send_design(f);
	return design_cohort_count(s);
}

size_t recv_model_count(const struct recv_model *model)
{
	assert(model);
	const struct design *design = recv_model_design(model);
	return design_count(design);
}

size_t recv_model_dim(const struct recv_model *m)
{
	assert(m);
	const struct frame *f = recv_model_frame(m);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	
	size_t dim = (design_trait_dim(r) + design_tvar_dim(r)
		      + design2_trait_dim(d) + design2_tvar_dim(d));
	return dim;
}

double recv_model_logsumwt0(const struct recv_model *m, size_t c)
{
	assert(m);
	assert(c < recv_model_cohort_count(m));
	return m->cohort_models[c].log_W0 + m->cohort_models[c].max_eta0;
}

double *recv_model_logwts0(const struct recv_model *m, size_t c)
{
	assert(m);
	assert(c < recv_model_cohort_count(m));
	return ((struct recv_model *)m)->cohort_models[c].eta0;
}

double *recv_model_probs0(const struct recv_model *m, size_t c)
{
	assert(m);
	assert(c < recv_model_cohort_count(m));
	return ((struct recv_model *)m)->cohort_models[c].p0;
}

double recv_model_prob0(const struct recv_model *m, size_t c, size_t jrecv)
{
	assert(m);
	assert(c < recv_model_cohort_count(m));
	assert(jrecv < recv_model_count(m));

	const double *p0 = m->cohort_models[c].p0;
	return p0[jrecv];
}

double *recv_model_mean0(const struct recv_model *m, size_t c)
{
	assert(m);
	assert(c < recv_model_cohort_count(m));
	return ((struct recv_model *)m)->cohort_models[c].mean0;
}

double *recv_model_imat0(const struct recv_model *m, size_t c)
{
	assert(m);
	assert(c < recv_model_cohort_count(m));
	return ((struct recv_model *)m)->cohort_models[c].imat0;
}

void recv_model_set_coefs(struct recv_model *m, const struct recv_coefs *coefs)
{
	const struct frame *f = recv_model_frame(m);
	size_t nc = recv_model_cohort_count(m);

	if (coefs) {
		memcpy(m->coefs.all, coefs->all, m->coefs.dim * sizeof(*m->coefs.all));
	} else {
		memset(m->coefs.all, 0, m->coefs.dim * sizeof(*m->coefs.all));		
	}

	size_t ic;
	for (ic = 0; ic < nc; ic++) {
		cohort_set(&m->cohort_models[ic], ic, f, coefs);
	}

	struct recv_model_sender *senders = m->sender_models;
	size_t isend, nsend = recv_model_send_count(m);
	for (isend = 0; isend < nsend; isend++) {
		sender_set(m, &senders[isend], isend, f, coefs);
	}
}

struct recv_model_sender *recv_model_send(const struct recv_model *m,
					  size_t isend)
{
	assert(m);

	struct recv_model_sender *rm =
	    sender_raw((struct recv_model *)m, isend);
	return rm;
}

void recv_model_get_active(const struct recv_model *m, size_t isend,
			   size_t **jrecv, size_t *n)
{
	assert(m);
	assert(isend < recv_model_send_count(m));
	assert(jrecv);
	assert(n);

	const struct recv_model_sender *sm = recv_model_send(m, isend);
	*jrecv = sm->active.indx;
	*n = sm->active.nz;
}

double recv_model_logsumwt(const struct recv_model *m, size_t isend)
{
	assert(m);
	assert(isend < recv_model_send_count(m));

	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->log_W + rm->scale;
}

double recv_model_invgrow(const struct recv_model *m, size_t isend)
{
	assert(m);
	assert(isend < recv_model_send_count(m));

	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->gamma;
}

/*
 * log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double recv_model_logprob(const struct recv_model *m, size_t isend,
			  size_t jrecv)
{
	assert(m);
	assert(isend < recv_model_send_count(m));
	assert(jrecv < recv_model_count(m));

	size_t c = recv_model_cohort(m, isend);
	const struct recv_model_cohort *cm = &m->cohort_models[c];
	const struct recv_model_sender *send = recv_model_send(m, isend);
	double deta = 0.0;
	ptrdiff_t iactive = vpattern_find(&send->active, jrecv);

	if (iactive >= 0) {
		deta = send->deta[iactive];
	}

	double scale = send->scale;
	double log_W = send->log_W;
	double eta0 = cm->eta0[jrecv];
	double eta = eta0 + deta;
	double log_p = (eta - scale) - log_W;

	log_p = MIN(0.0, log_p);

	assert(log_p <= 0.0);
	return log_p;
}

double recv_model_prob(const struct recv_model *m, size_t isend, size_t jrecv)
{
	assert(m);
	assert(isend < recv_model_send_count(m));
	assert(jrecv < recv_model_count(m));

	double lp = recv_model_logprob(m, isend, jrecv);
	double p = exp(lp);
	return p;
}

void recv_model_axpy_probs(double alpha,
			   const struct recv_model *m,
			   size_t isend, double *y)
{
	assert(m);
	assert(isend < recv_model_send_count(m));

	size_t j, n = recv_model_count(m);

	for (j = 0; j < n; j++) {
		double p = recv_model_prob(m, isend, j);
		y[j] += alpha * p;
	}
}

const size_t *recv_model_cohorts(const struct recv_model *model)
{
	const struct frame *f = recv_model_frame(model);
	const struct design *s = frame_send_design(f);
	const size_t *cs, *reps;
	size_t nc;
	
	design_get_cohorts(s, &cs, &reps, &nc);
	return cs;
}
