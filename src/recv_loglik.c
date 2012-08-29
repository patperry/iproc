#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "xalloc.h"
#include "recv_loglik.h"

static void cohort_init(struct recv_loglik_cohort *cll);
static void cohort_deinit(struct recv_loglik_cohort *cll);
static void cohort_clear(struct recv_loglik_cohort *cll);
static void cohort_add(struct recv_loglik_cohort *cll);

static void sender_init(struct recv_loglik_sender *sll, const struct frame *f);
static void sender_deinit(struct recv_loglik_sender *sll);
static void sender_add(struct recv_loglik_sender *sll, const struct frame *f,
		       const struct mlogitaug *m1, size_t isend,
		       const size_t *jrecv, size_t nto,
		       struct recv_loglik_update *last);
static void sender_clear(struct recv_loglik_sender *sll);


void cohort_init(struct recv_loglik_cohort *cll)
{
}

void cohort_deinit(struct recv_loglik_cohort *cll)
{

}

void cohort_clear(struct recv_loglik_cohort *cll)
{
}

void cohort_add(struct recv_loglik_cohort *cll)
{
}



void sender_init(struct recv_loglik_sender *sll, const struct frame *f)
{
	recv_coefs_init(&sll->mean, f);
	sender_clear(sll);
}

void sender_deinit(struct recv_loglik_sender *sll)
{
	recv_coefs_deinit(&sll->mean);
}

void sender_clear(struct recv_loglik_sender *sll)
{
	size_t dim = sll->mean.dim;
	sll->count = 0;
	sll->dev = 0;
	memset(sll->mean.all, 0, dim * sizeof(*sll->mean.all));
}


void sender_add(struct recv_loglik_sender *sll, const struct frame *f,
		const struct mlogitaug *m1, size_t isend,
		const size_t *jrecv, size_t nto,
		struct recv_loglik_update *last)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);

	size_t base_dim = mlogit_dim(mlogitaug_base(m1));
	size_t dim1 = mlogitaug_dim(m1);
	size_t dim = base_dim + dim1;

	const struct catdist1 *dist = mlogitaug_cached_dist(m1);

	double dev = 0.0;
	size_t count = nto;
	size_t ito;

	/* compute dev and set score := observed */
	memset(last->score.all, 0, dim * sizeof(*last->score.all));
	for (ito = 0; ito < nto; ito++) {
		double lp = catdist1_cached_lprob(dist, jrecv[ito]);
		dev += -2 * lp;
		design_axpy(1.0, r, jrecv[ito], &last->score.recv);
		design2_axpy(1.0, d, isend, jrecv[ito], &last->score.dyad);
	}

	/* set count, dev */
	last->count = count;
	last->dev = dev;

	/* compute mean */
	const double *base_mean = mlogitaug_cached_base_mean(m1);
	const double *mean = mlogitaug_cached_mean(m1);
	blas_dcopy(base_dim, base_mean, 1, last->mean.all, 1);
	blas_dcopy(dim1, mean, 1, last->mean.all + base_dim, 1);

	/* compute score := observed - expected */
	blas_daxpy(dim, -(double)nto, last->mean.all, 1, last->score.all, 1);

	/* update totals */
	sll->count += last->count;
	sll->dev += last->dev;
	blas_daxpy(dim, last->count, last->mean.all, 1, sll->mean.all, 1);
}


void recv_loglik_init(struct recv_loglik *ll, struct recv_model *m)
{
	assert(ll);
	assert(m);

	const struct frame *f = recv_model_frame(m);

	ll->model = m;

	size_t ic, nc = recv_model_cohort_count(m);
	struct recv_loglik_cohort *cohorts = xcalloc(nc, sizeof(*cohorts));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cohorts[ic]);
	}
	ll->cohorts = cohorts;

	size_t is, ns = recv_model_send_count(m);
	struct recv_loglik_sender *senders = xcalloc(ns, sizeof(*senders));
	for (is = 0; is < ns; is++) {
		sender_init(&senders[is], f);
	}
	ll->senders = senders;

	recv_coefs_init(&ll->last.mean, recv_model_frame(m));
	recv_coefs_init(&ll->last.score, recv_model_frame(m));

	recv_loglik_clear(ll);
}


void recv_loglik_deinit(struct recv_loglik *ll)
{
	assert(ll);

	recv_coefs_deinit(&ll->last.score);
	recv_coefs_deinit(&ll->last.mean);

	struct recv_loglik_sender *senders = ll->senders;
	size_t is, ns = recv_model_send_count(ll->model);
	for (is = 0; is < ns; is++) {
		sender_deinit(&senders[is]);
	}
	free(senders);

	struct recv_loglik_cohort *cohorts = ll->cohorts;
	size_t ic, nc = recv_model_cohort_count(ll->model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cohorts[ic]);
	}
	free(cohorts);
}


void recv_loglik_clear(struct recv_loglik *ll)
{
	assert(ll);

	const struct recv_model *m = ll->model;
	struct recv_loglik_cohort *cohorts = ll->cohorts;
	size_t ic, nc = recv_model_cohort_count(m);

	for (ic = 0; ic < nc; ic++) {
		cohort_clear(&cohorts[ic]);
	}

	struct recv_loglik_sender *senders = ll->senders;
	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sender_clear(&senders[is]);
	}

	ll->last.count = 0;
	ll->last.dev = 0;
	memset(ll->last.mean.all, 0, ll->last.mean.dim * sizeof(*ll->last.mean.all));
	memset(ll->last.score.all, 0, ll->last.score.dim * sizeof(*ll->last.score.all));
}


void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f, const struct message *msg)
{
	size_t isend = msg->from;
	const struct mlogitaug *m1 = recv_model_mlogit(ll->model, isend);
	size_t c = recv_model_cohort(ll->model, isend);
	sender_add(&ll->senders[isend], f, m1, msg->from, msg->to, msg->nto, &ll->last);
	cohort_add(&ll->cohorts[c]);
}


void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f, const struct messages *msgs)
{
	struct messages_iter it;
	const struct message *msg;
	size_t i, n;

	MESSAGES_FOREACH(it, msgs) {
		double t = MESSAGES_TIME(it);
		frame_advance(f, t);

		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			recv_loglik_add(ll, f, msg);
		}

		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(f, msg);
		}
	}

}


size_t recv_loglik_count(const struct recv_loglik *ll)
{
	size_t count = 0;
	const struct recv_model *m = ll->model;
	const struct recv_loglik_sender *sll;

	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sll = &ll->senders[is];
		count += sll->count;
	}

	return count;
}


double recv_loglik_dev(const struct recv_loglik *ll)
{
	double dev = 0;
	const struct recv_model *m = ll->model;
	const struct recv_loglik_sender *sll;

	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sll = &ll->senders[is];
		dev += sll->dev;
	}

	return dev;
}


void recv_loglik_axpy_mean(double alpha, const struct recv_loglik *ll, struct recv_coefs *y)
{
	const struct recv_model *m = ll->model;
	const struct recv_loglik_sender *sll;
	size_t dim = recv_model_dim(m);

	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sll = &ll->senders[is];
		blas_daxpy(dim, alpha, sll->mean.all, 1, y->all, 1);
	}
}

size_t recv_loglik_last_count(const struct recv_loglik *ll)
{
	return ll->last.count;
}


double recv_loglik_last_dev(const struct recv_loglik *ll)
{
	return ll->last.dev;
}

void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll, struct recv_coefs *y)
{
	assert(ll->last.mean.dim == y->dim);
	double scale = ll->last.count * alpha;
	blas_daxpy(ll->last.mean.dim, scale, ll->last.mean.all, 1, y->all, 1);
}

void recv_loglik_axpy_last_score(double alpha, const struct recv_loglik *ll, struct recv_coefs *y)
{
	assert(ll->last.score.dim == y->dim);
	blas_daxpy(ll->last.score.dim, alpha, ll->last.score.all, 1, y->all, 1);
}
