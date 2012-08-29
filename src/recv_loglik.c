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

static void sender_init(struct recv_loglik_sender *sll);
static void sender_deinit(struct recv_loglik_sender *sll);
void sender_add(struct recv_loglik_sender *sll, const struct mlogitaug *m1,
		const size_t *jrecv, size_t nto, struct recv_loglik_update *last);
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



void sender_init(struct recv_loglik_sender *sll)
{
	sender_clear(sll);
}

void sender_deinit(struct recv_loglik_sender *sll)
{
}

void sender_clear(struct recv_loglik_sender *sll)
{
	sll->count = 0;
	sll->dev = 0;
}


void sender_add(struct recv_loglik_sender *sll, const struct mlogitaug *m1,
		const size_t *jrecv, size_t nto, struct recv_loglik_update *last)
{
	size_t base_dim = mlogit_dim(mlogitaug_base(m1));
	size_t dim = mlogitaug_dim(m1);

	const struct catdist1 *dist = mlogitaug_cached_dist(m1);
	const double *base_mean = mlogitaug_cached_base_mean(m1);
	const double *mean = mlogitaug_cached_mean(m1);

	double dev = 0.0;
	size_t count = nto;
	size_t ito;

	for (ito = 0; ito < nto; ito++) {
		double lp = catdist1_cached_lprob(dist, jrecv[ito]);
		dev += -2 * lp;
	}

	last->count = count;
	last->dev = dev;

	blas_dcopy(base_dim, base_mean, 1, last->mean.all, 1);
	blas_dcopy(dim, mean, 1, last->mean.all + base_dim, 1);

	sll->count += last->count;
	sll->dev += last->dev;
}


void recv_loglik_init(struct recv_loglik *ll, struct recv_model *m)
{
	assert(ll);
	assert(m);

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
		sender_init(&senders[is]);
	}
	ll->senders = senders;

	recv_coefs_init(&ll->last.mean, recv_model_frame(m));

	recv_loglik_clear(ll);
}


void recv_loglik_deinit(struct recv_loglik *ll)
{
	assert(ll);

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
}


void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f, const struct message *msg)
{
	size_t isend = msg->from;
	const struct mlogitaug *m1 = recv_model_mlogit(ll->model, isend);
	size_t c = recv_model_cohort(ll->model, isend);
	sender_add(&ll->senders[isend], m1, msg->to, msg->nto, &ll->last);
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

