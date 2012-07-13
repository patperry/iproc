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
#include "sblas.h"
#include "util.h"
#include "recv_model.h"

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
		sm->W = exp(sm->log_W);
	}

out:
	sm->gamma = exp((log_W0 + (max_eta0 - sm->scale)) - sm->log_W);
	assert(sm->gamma >= 0.0);
	assert(isfinite(sm->log_W));
}

static void cohort_set(struct recv_model_cohort *cm, size_t isend,
		       const struct recv_frame *rf,
		       const struct recv_coefs *coefs)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	size_t dim = recv_fmla_trait_dim(fmla);

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
	recv_frame_mul0(1.0, rf, isend, coefs->traits, 0.0, cm->eta0);

	/* max_eta0 */
	cm->max_eta0 = vector_max(nrecv, cm->eta0);

	/* store log_p0 in p0 */
	blas_dcopy(nrecv, cm->eta0, 1, cm->p0, 1);
	vector_shift(nrecv, -cm->max_eta0, cm->p0);	/* guard against overflow */

	/* log_W0 */
	cm->log_W0 = vector_logsumexp(nrecv, cm->p0);
	vector_shift(nrecv, -cm->log_W0, cm->p0);
	vector_exp(nrecv, cm->p0);

	/* mean0 */
	recv_frame_tmul0(1.0, rf, isend, cm->p0, 0.0, cm->mean0);

	/* imat0 */
	double pj;
	size_t jrecv;

	double *y = xmalloc(dim * sizeof(double));
	matrix_dzero(dim, dim, &cm->imat0);

	for (jrecv = 0; jrecv < nrecv; jrecv++) {
		blas_dcopy(dim, cm->mean0, 1, y, 1);
		recv_frame_axpy0(-1.0, rf, isend, jrecv, y);
		pj = cm->p0[jrecv];

		blas_dger(dim, dim, pj, y, 1, y, 1, &cm->imat0);
	}
	free(y);

#ifndef NDEBUG
	/* W0 */
	cm->W0 = exp(cm->log_W0 + cm->max_eta0);

	/* w0 */
	blas_dcopy(nrecv, cm->p0, 1, cm->w0, 1);
	blas_dscal(nrecv, cm->W0, cm->w0, 1);
#endif
}

static void cohort_init(struct recv_model_cohort *cm,
			size_t isend,			
			const struct recv_frame *rf,
			const struct recv_coefs *coefs)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	size_t dim = recv_fmla_trait_dim(fmla);
	size_t dim2 = dim * dim;

	cm->eta0 = xmalloc(nrecv * sizeof(*cm->eta0));
	cm->p0 = xmalloc(nrecv * sizeof(*cm->p0));
	cm->mean0 = xmalloc(dim * sizeof(*cm->mean0));
	cm->imat0.data = xmalloc(dim2 * sizeof(double));
	cm->imat0.lda = MAX(1, dim);
#ifndef NDEBUG
	cm->w0 = xmalloc(nrecv * sizeof(*cm->w0));
#endif
	cohort_set(cm, isend, rf, coefs);
}

static void cohort_deinit(struct recv_model_cohort *cm)
{
	assert(cm);
#ifndef NDEBUG
	free(cm->w0);
#endif
	free(cm->imat0.data);
	free(cm->mean0);
	free(cm->p0);
	free(cm->eta0);
}


static void sender_clear(const struct recv_model *m,
			 struct recv_model_sender *send, size_t isend)
{
	const struct recv_fmla *fmla = recv_model_fmla(m);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	size_t c = recv_fmla_cohort(fmla, isend);
	const struct recv_model_cohort *cm = &m->cohorts[c];
	double max_eta0 = cm->max_eta0;
	double log_W0 = cm->log_W0;
	double W0 = cm->W0;

	vpattern_clear(&send->active);

	/* take the initial values if there are self-loops */
	if (frame_has_loops(f)) {
		send->gamma = 1.0;
		send->scale = max_eta0;
		send->log_W = log_W0;
		send->W = W0;
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
		       const struct recv_coefs *coefs)
{
	const struct recv_frame *rf = &m->recv_frame;
	const struct recv_fmla *fmla = recv_model_fmla(m);
	const struct frame *f = recv_fmla_frame(fmla);
	double *zbuf = m->zbuf;
	
	sender_clear(m, send, isend);

	int has_loops = frame_has_loops(f);

	/* compute the eta values */
	const size_t *active = send->active.indx;
	size_t iz, nz = send->active.nz;

	recv_frame_mul1(1.0, rf, isend, coefs->tvars, 1.0, zbuf);
	for (iz = 0; iz < nz; iz++) {
		size_t ix = active[iz];
		send->deta[iz] = zbuf[ix];
		zbuf[ix] = 0;
	}

	if (!has_loops) {
		ptrdiff_t ix = vpattern_find(&send->active, isend);
		assert(ix >= 0);
		send->deta[ix] = -INFINITY;
	}

	/* compute the changes in weights */
	size_t c = recv_fmla_cohort(fmla, isend);
	const struct recv_model_cohort *cm = &m->cohorts[c];
	size_t nrecv = frame_recv_count(f);
	compute_weight_changes(send, cm, nrecv);

	assert(send->gamma >= 0.0);
	assert(isfinite(send->log_W));
}

static void sender_init(const struct recv_model *m,
			struct recv_model_sender *send, size_t isend)
{
	assert(send);
	assert(m);

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

	return &m->senders[isend];
}

static void recv_model_recv_update(void *udata, struct recv_frame *rf, size_t isend,
				   size_t jrecv, const double *delta,
				   const struct vpattern *pat)
{
	struct recv_model *m = udata;
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	
	if (isend == jrecv && !frame_has_loops(f))
		return;

	size_t icohort = recv_fmla_cohort(fmla, isend);
	struct recv_model_sender *send = sender_raw(m, isend);
	const struct recv_model_cohort *cm = &m->cohorts[icohort];
	const struct recv_coefs *coefs = &m->coefs;

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

	const double *dyn_coefs = coefs->tvars;

	double deta = sblas_ddoti(delta, pat, dyn_coefs);
	double scale = send->scale;
	double gamma = send->gamma;
	double log_W = send->log_W;
	double W = send->W;
	double eta = cm->eta0[jrecv] + (*pdeta);
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
		sender_set(m, send, isend, coefs);
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
	const struct recv_frame *rf = &m->recv_frame;
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	struct recv_model_sender *senders = m->senders;
	size_t i, n = frame_send_count(f);

	for (i = 0; i < n; i++) {
		sender_clear(m, &senders[i], i);
	}
}

static void recv_model_clear(void *udata, struct recv_frame *rf)
{
	(void)rf;		// unused
	struct recv_model *m = udata;
	assert(m);
	model_clear(m);
}

void recv_model_init(struct recv_model *model,
		     const struct recv_fmla *fmla,
		     const struct recv_coefs *coefs)
{
	assert(coefs->fmla == fmla);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nsend = frame_send_count(f);
	
	recv_frame_init(&model->recv_frame, fmla);
	recv_coefs_init_copy(&model->coefs, coefs);
	model->zbuf = xcalloc(nsend, sizeof(*model->zbuf));
	
	const size_t *cs, *ireps;
	size_t ic, nc;
	recv_fmla_get_cohorts(fmla, &cs, &ireps, &nc);
	
	struct recv_model_cohort *cohorts = xcalloc(nc, sizeof(*cohorts));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cohorts[ic], ireps[ic], &model->recv_frame, &model->coefs);
	}
	model->cohorts = cohorts;

	size_t isend;
	struct recv_model_sender *sms = xcalloc(nsend, sizeof(*sms));

	for (isend = 0; isend < nsend; isend++) {
		sender_init(model, &sms[isend], isend);
	}
	model->senders = sms;

	struct recv_frame_callbacks callbacks = {
		recv_model_recv_update,
		NULL,
		recv_model_clear
	};

	recv_frame_add_observer(&model->recv_frame, model, &callbacks);
}

void recv_model_deinit(struct recv_model *model)
{
	const struct recv_frame *rf = &model->recv_frame;
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	
	recv_frame_remove_observer(&model->recv_frame, model);

	struct recv_model_sender *sms = model->senders;
	size_t isend, nsend = frame_send_count(f);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&sms[isend]);
	}
	free(sms);

	struct recv_model_cohort *cms = model->cohorts;
	size_t ic, nc = recv_fmla_cohort_count(fmla);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cms[ic]);
	}
	free(cms);

	free(model->zbuf);
	recv_coefs_deinit(&model->coefs);
	recv_frame_deinit(&model->recv_frame);
}


double recv_model_logsumwt0(const struct recv_model *m, size_t c)
{
	assert(m);
	return m->cohorts[c].log_W0 + m->cohorts[c].max_eta0;
}

double *recv_model_logwts0(const struct recv_model *m, size_t c)
{
	assert(m);
	return ((struct recv_model *)m)->cohorts[c].eta0;
}

double *recv_model_probs0(const struct recv_model *m, size_t c)
{
	assert(m);
	return ((struct recv_model *)m)->cohorts[c].p0;
}

double recv_model_prob0(const struct recv_model *m, size_t c, size_t jrecv)
{
	assert(m);

	const double *p0 = m->cohorts[c].p0;
	return p0[jrecv];
}

double *recv_model_mean0(const struct recv_model *m, size_t c)
{
	assert(m);
	return ((struct recv_model *)m)->cohorts[c].mean0;
}

struct dmatrix *recv_model_imat0(const struct recv_model *m, size_t c)
{
	assert(m);
	return &((struct recv_model *)m)->cohorts[c].imat0;
}

void recv_model_set_coefs(struct recv_model *m, const struct recv_coefs *coefs)
{
	const struct recv_frame *rf = &m->recv_frame;
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nc = recv_fmla_cohort_count(fmla);

	if (coefs) {
		recv_coefs_assign_copy(&m->coefs, coefs);
	} else {
		recv_coefs_clear(&m->coefs);
	}

	size_t ic;
	for (ic = 0; ic < nc; ic++) {
		size_t isend = recv_fmla_cohort_rep(fmla, ic);
		cohort_set(&m->cohorts[ic], isend, rf, &m->coefs);
	}

	struct recv_model_sender *senders = m->senders;
	size_t isend, nsend = frame_send_count(f);
	for (isend = 0; isend < nsend; isend++) {
		sender_set(m, &senders[isend], isend, &m->coefs);
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
	const struct recv_model_sender *sm = recv_model_send(m, isend);
	*jrecv = sm->active.indx;
	*n = sm->active.nz;
}

double recv_model_logsumwt(const struct recv_model *m, size_t isend)
{
	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->log_W + rm->scale;
}

double recv_model_invgrow(const struct recv_model *m, size_t isend)
{
	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->gamma;
}

/*
 * log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double recv_model_logprob(const struct recv_model *m, size_t isend,
			  size_t jrecv)
{
	const struct recv_frame *rf = &m->recv_frame;
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	size_t c = recv_fmla_cohort(fmla, isend);
	const struct recv_model_cohort *cm = &m->cohorts[c];
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
	double lp = recv_model_logprob(m, isend, jrecv);
	double p = exp(lp);
	return p;
}

void recv_model_axpy_probs(double alpha,
			   const struct recv_model *m,
			   size_t isend, double *y)
{
	const struct recv_frame *rf = &m->recv_frame;
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t j, n = frame_recv_count(f);

	for (j = 0; j < n; j++) {
		double p = recv_model_prob(m, isend, j);
		y[j] += alpha * p;
	}
}

