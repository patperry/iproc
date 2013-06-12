#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "xalloc.h"
#include "recv_loglik.h"

static void sender_init(struct recv_loglik_sender *sll, const struct recv_model *m);
static void sender_deinit(struct recv_loglik_sender *sll);
static void sender_add(struct recv_loglik_sender *sll, const struct design *r,
		       const struct design2 *d, const struct mlogitaug *m1,
		       size_t isend, const size_t *jrecv, size_t nto,
		       struct recv_loglik_update *last);
static void sender_clear(struct recv_loglik_sender *sll, const struct design *r, const struct design2 *d);

#ifndef NDEBUG
static int recv_loglik_moments(const struct recv_loglik *ll);
#endif


void sender_init(struct recv_loglik_sender *sll, const struct recv_model *m)
{
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	
	recv_params_init(&sll->mean, r, d);
	recv_params_init(&sll->score, r, d);

	size_t dim = recv_model_dim(m);
	size_t cov_dim = dim * (dim + 1) / 2;
	sll->cov = xmalloc(cov_dim * sizeof(*sll->cov));

	sender_clear(sll, r, d);
}

void sender_deinit(struct recv_loglik_sender *sll)
{
	free(sll->cov);
	recv_params_deinit(&sll->score);
	recv_params_deinit(&sll->mean);
}

void sender_clear(struct recv_loglik_sender *sll, const struct design *r, const struct design2 *d)
{
	size_t dim = design_dim(r) + design2_dim(d);
	size_t cov_dim = dim * (dim + 1) / 2;
	sll->count = 0;
	sll->dev = 0;
	recv_params_set(&sll->mean, NULL, r, d);
	recv_params_set(&sll->score, NULL, r, d);
	memset(sll->cov, 0, cov_dim * sizeof(*sll->cov));
}


static void copy_cov(const struct mlogitaug *m1, double *dst)
{
	size_t base_dim = mlogit_dim(mlogitaug_base(m1));
	size_t aug_dim = mlogitaug_dim(m1);
	//size_t dim = base_dim + aug_dim;
	size_t base_cov_dim = base_dim * (base_dim + 1) / 2;
	size_t aug_cov_dim = aug_dim * (aug_dim + 1) / 2;
	const double *base_cov = mlogitaug_base_cov(m1);
	const double *cross_cov = mlogitaug_cross_cov(m1);
	const double *aug_cov = mlogitaug_cov(m1);
	size_t i;

	assert(RECV_LOGLIK_UPLO == MLOGIT_COV_UPLO);
	if (MLOGIT_COV_UPLO == BLAS_UPPER) {
		for (i = 0; i < base_dim; i++) {
			size_t rowlen = base_dim - i;
			memcpy(dst, base_cov, rowlen * sizeof(double));
			dst += rowlen;
			base_cov += rowlen;

			memcpy(dst, cross_cov, aug_dim * sizeof(double));
			dst += aug_dim;
			cross_cov += aug_dim;
		}
		memcpy(dst, aug_cov, aug_cov_dim * sizeof(double));
	} else {
		memcpy(dst, base_cov, base_cov_dim * sizeof(double));
		dst += base_cov_dim;

		for (i = 0; i < aug_dim; i++) {
			memcpy(dst, cross_cov, base_dim * sizeof(double));
			dst += base_dim;
			cross_cov += base_dim;

			size_t rowlen = i + 1;
			memcpy(dst, aug_cov, rowlen * sizeof(double));
			dst += rowlen;
			aug_cov += rowlen;
		}
	}
}

void sender_add(struct recv_loglik_sender *sll, const struct design *r,
		const struct design2 *d,
		const struct mlogitaug *m1, size_t isend,
		const size_t *jrecv, size_t nto,
		struct recv_loglik_update *last)
{
	int moments = mlogit_moments(mlogitaug_base(m1));

	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_tvar_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dimd1 = design2_tvar_dim(d);

	size_t base_dim = mlogit_dim(mlogitaug_base(m1));
	size_t aug_dim = mlogitaug_dim(m1);
	size_t dim = base_dim + aug_dim;
	size_t cov_dim = dim * (dim + 1) / 2;

	const struct catdist1 *dist = mlogitaug_dist(m1);

	double dev = 0.0;
	size_t count = nto;
	size_t ito;

	/* compute dev and set score := observed */
	if (moments >= 1) {
		recv_params_set(&last->score, NULL, r, d);
	}
	for (ito = 0; ito < nto; ito++) {
		double lp = catdist1_lprob(dist, jrecv[ito]);
		dev += -2 * lp;
		if (moments >= 1) {
			design_axpy(1.0, r, jrecv[ito], &last->score.recv);
			design2_axpy(1.0, d, isend, jrecv[ito], &last->score.dyad);
		}
	}

	/* set count, dev */
	last->count = count;
	last->dev = dev;

	sll->count += last->count;
	sll->dev += last->dev;

	if (moments < 1)
		return;

	/* compute mean */
	const double *base_mean = mlogitaug_base_mean(m1);
	const double *mean = mlogitaug_mean(m1);

	//blas_dcopy(base_dim, base_mean, 1, last->mean.all, 1);
	memcpy(last->mean.recv.traits, base_mean, dimr0 * sizeof(double));
	memcpy(last->mean.recv.tvars, base_mean + dimr0, dimr1 * sizeof(double));
	memcpy(last->mean.dyad.traits, base_mean + dimr0 + dimr1, dimd0 * sizeof(double));
	
	//blas_dcopy(aug_dim, mean, 1, last->mean.all + base_dim, 1);
	memcpy(last->mean.dyad.tvars, mean, dimd1 * sizeof(double));


	/* compute score := observed - expected */
	recv_params_axpy(-(double)nto, &last->mean, &last->score, r, d);

	recv_params_axpy(last->count, &last->mean, &sll->mean, r, d);
	recv_params_axpy(1.0, &last->score, &sll->score, r, d);

	if (moments < 2)
		return;

	/* copy cov */
	copy_cov(m1, last->cov);

	blas_daxpy(cov_dim, last->count, last->cov, 1, sll->cov, 1);
}


void recv_loglik_init(struct recv_loglik *ll, struct recv_model *m)
{
	assert(ll);
	assert(m);

	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);

	ll->model = m;

	size_t is, ns = recv_model_send_count(m);
	struct recv_loglik_sender *senders = xcalloc(ns, sizeof(*senders));
	for (is = 0; is < ns; is++) {
		sender_init(&senders[is], m);
	}
	ll->senders = senders;

	recv_params_init(&ll->last.mean, r, d);
	recv_params_init(&ll->last.score, r, d);

	size_t dim = recv_model_dim(m);
	size_t cov_dim = dim * (dim + 1) / 2;
	ll->last.cov = xmalloc(cov_dim * sizeof(*ll->last.cov));

	recv_loglik_clear(ll);
}


void recv_loglik_deinit(struct recv_loglik *ll)
{
	assert(ll);

	free(ll->last.cov);
	recv_params_deinit(&ll->last.score);
	recv_params_deinit(&ll->last.mean);

	struct recv_loglik_sender *senders = ll->senders;
	size_t is, ns = recv_model_send_count(ll->model);
	for (is = 0; is < ns; is++) {
		sender_deinit(&senders[is]);
	}
	free(senders);
}


void recv_loglik_clear(struct recv_loglik *ll)
{
	assert(ll);

	const struct recv_model *m = ll->model;
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	struct recv_loglik_sender *senders = ll->senders;
	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sender_clear(&senders[is], r, d);
	}

	ll->last.count = 0;
	ll->last.dev = 0;
	recv_params_set(&ll->last.mean, NULL, r, d);
	recv_params_set(&ll->last.score, NULL, r, d);
}


void recv_loglik_add(struct recv_loglik *ll, const struct message *msg)
{
	const struct recv_model *m = recv_loglik_model(ll);
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	struct history *hr = design_history(r);
	struct history *hd = design2_history(d);
	double t = msg->time;
	size_t isend = msg->from;

	if (t < history_time(hr)) {
		history_reset(hr);
	}
	history_advance(hr, t);

	if (t < history_time(hd)) {
		history_reset(hd);
	}
	history_advance(hd, t);

	const struct mlogitaug *m1 = recv_model_mlogit(ll->model, isend);
	sender_add(&ll->senders[isend], r, d, m1, msg->from, msg->to, msg->nto, &ll->last);
}


void recv_loglik_add_all(struct recv_loglik *ll, const struct message *msgs, size_t n)
{
	size_t i;

	for (i = 0; i < n; i++) {
		const struct message *msg = &msgs[i];
		recv_loglik_add(ll, msg);
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


void recv_loglik_axpy_mean(double alpha, const struct recv_loglik *ll, struct recv_params *y)
{
	assert(recv_loglik_moments(ll) >= 1);
	const struct recv_model *m = ll->model;
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	const struct recv_loglik_sender *sll;

	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sll = &ll->senders[is];
		recv_params_axpy(alpha, &sll->mean, y, r, d);
	}
}


void recv_loglik_axpy_score(double alpha, const struct recv_loglik *ll, struct recv_params *y)
{
	assert(recv_loglik_moments(ll) >= 1);
	const struct recv_model *m = ll->model;
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	const struct recv_loglik_sender *sll;

	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sll = &ll->senders[is];
		recv_params_axpy(alpha, &sll->score, y, r, d);
	}
}


void recv_loglik_axpy_imat(double alpha, const struct recv_loglik *ll, double *y)
{
	assert(recv_loglik_moments(ll) >= 2);
	const struct recv_model *m = ll->model;
	const struct recv_loglik_sender *sll;
	size_t dim = recv_model_dim(m);
	size_t cov_dim = dim * (dim + 1) / 2;

	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sll = &ll->senders[is];
		blas_daxpy(cov_dim, alpha, sll->cov, 1, y, 1);
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

void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll, struct recv_params *y)
{
	assert(recv_loglik_moments(ll) >= 1);

	const struct recv_model *m = recv_loglik_model(ll);
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	double scale = ll->last.count * alpha;

	recv_params_axpy(scale, &ll->last.mean, y, r, d);
}

void recv_loglik_axpy_last_score(double alpha, const struct recv_loglik *ll, struct recv_params *y)
{
	assert(recv_loglik_moments(ll) >= 1);

	const struct recv_model *m = recv_loglik_model(ll);
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);

	recv_params_axpy(alpha, &ll->last.score, y, r, d);
}

void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll, double *y)
{
	assert(recv_loglik_moments(ll) >= 2);
	size_t dim = recv_model_dim(ll->model);
	size_t cov_dim = dim * (dim + 1) / 2;
	double scale = ll->last.count * alpha;
	blas_daxpy(cov_dim, scale, ll->last.cov, 1, y, 1);
}

#ifndef NDEBUG
int recv_loglik_moments(const struct recv_loglik *ll)
{
	return recv_model_moments(ll->model);
}
#endif
