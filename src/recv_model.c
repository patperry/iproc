#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "lapack.h"
#include "matrixutil.h"
#include "util.h"
#include "sblas.h"
#include "recv_model.h"


static size_t cohort_count(const struct recv_model *m)
{
	const struct design2 *d = recv_model_design2(m);
	return design2_cohort_count(d);
}


static void cohort_clear(struct recv_model_cohort *cm, size_t c,
			 const struct design *r, const struct design2 *d)
{
	const struct history *h = design_history(r);
	size_t isend = design2_cohort_rep(d, c);
	size_t nrecv = design_count(r);
	size_t dimr = design_dim(r);
	size_t dimr0 = design_trait_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dim = dimr + dimd0;
	const double *xr = design_trait_matrix(r);
	const double *xd = design2_trait_matrix(d, isend);
	double *x = xcalloc(nrecv * dim, sizeof(*x));

	if (dimr0 > 0) {
		assert(dim > 0);
		lapack_dlacpy(LA_COPY_ALL, dimr0, nrecv, xr, dimr0, x, dim);
	}

	if (dimd0 > 0) {
		assert(dim > 0);
		lapack_dlacpy(LA_COPY_ALL, dimd0, nrecv, xd, dimd0, x + dimr, dim);
	}

	mlogit_set_all_x(&cm->mlogit, x);

	version_watch_set(&cm->version, history_version(h));
	cm->tcur = -INFINITY;

	free(x);
}


static void cohort_set(struct recv_model_cohort *cm, size_t c,
		       const struct design *r, const struct design2 *d,
		       const struct recv_params *p)
{
	size_t dimr = design_dim(r);
	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_trait_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dim = dimr + dimd0;
	double *beta = xcalloc(dim, sizeof(double));

	if (p) {
		if (p->recv.traits) {
			memcpy(beta, p->recv.traits, dimr0 * sizeof(double));
		}
		if (p->recv.tvars) {
			memcpy(beta + dimr0, p->recv.tvars, dimr1 * sizeof(double));
		}
		if (p->dyad.traits) {
			memcpy(beta + dimr, p->dyad.traits, dimd0 * sizeof(double));
		}
	}

	mlogit_set_coefs(&cm->mlogit, beta);

	free(beta);
}


static void cohort_update(struct recv_model_cohort *cm, size_t c,
			  const struct design *r, const struct design2 *d)
{
	const struct history *h = design_history(r);
	const struct version *v = history_version(h);
	const struct deltaset *ds = design_changes(r);
	const struct delta *delta;
	size_t j;
	const double *xj;
	const double *x;
	const size_t *ind;
	size_t off = design_trait_dim(r);
	size_t dim = design_tvar_dim(r);
	size_t nz;
	ptrdiff_t iz;

	if (version_changed(v, &cm->version)) {
		cohort_clear(cm, c, r, d);
	}

	double t0 = cm->tcur;
	double t = history_time(h);

	if (t0 == t)
		return;

	design_get_tvar_matrix(r, &x, &ind, &nz);

	for (delta = deltaset_head(ds); delta != NULL; delta = delta->prev) {
		if (delta->t < t0 || delta->t == -INFINITY)
			break;

		j = delta->i;
		iz = find_index(j, ind, nz);
		assert(iz >= 0);
		xj = x + iz * dim;

		mlogit_set_x(&cm->mlogit, j, off, dim, xj);
	}

	cm->tcur = t;
}


static void cohort_set_moments(struct recv_model_cohort *cm, int k)
{
	mlogit_set_moments(&cm->mlogit, k);
}


static void cohort_init(struct recv_model_cohort *cm,
			size_t c,
			const struct design *r, const struct design2 *d,
			const struct recv_params *p,
			struct mlogit_work *work)
{
	struct history *h = design_history(r);
	size_t nrecv = design_count(r);
	size_t dimr = design_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dim = dimr + dimd0;

	assert(design2_history(d) == h);

	mlogit_init(&cm->mlogit, nrecv, dim, work);
	version_watch_init(&cm->version, history_version(h));
	cohort_clear(cm, c, r, d);
	cohort_set(cm, c, r, d, p);
}


static void cohort_deinit(struct recv_model_cohort *cm, struct history *h)
{
	version_watch_deinit(&cm->version, history_version(h));
	mlogit_deinit(&cm->mlogit);
}


size_t recv_model_cohort(const struct recv_model *m, size_t isend)
{
	const struct design2 *d = recv_model_design2(m);
	return design2_cohort(d, isend);
}


static void sender_clear(struct recv_model_sender *sm, const struct design2 *d)
{
	const struct history *h = design2_history(d);
	mlogitaug_set_all_x(&sm->mlogitaug, NULL, NULL, 0);
	version_watch_set(&sm->version, history_version(h));
	sm->tcur = -INFINITY;
}


static void sender_set(struct recv_model_sender *sm,
		       size_t isend,
		       const struct design2 *d, const struct recv_params *p)
{
	const double *coefs = p ? p->dyad.tvars : NULL;
	mlogitaug_set_coefs(&sm->mlogitaug, coefs);
	//mlogitaug_check(&sm->mlogitaug);
}


static void sender_set_exclude_loops(struct recv_model_sender *sm, size_t isend,
				     int exclude_loops)
{
	if (exclude_loops) {
		const double neginf = -INFINITY;
		mlogitaug_set_all_offset(&sm->mlogitaug, &isend, &neginf, 1);
	} else {
		mlogitaug_set_all_offset(&sm->mlogitaug, NULL, NULL, 0);
	}
}


static void sender_update(struct recv_model_sender *sm, size_t isend,
			  const struct design2 *d)
{
	const struct history *h = design2_history(d);
	const struct version *v = history_version(h);
	const struct deltaset *ds = design2_changes(d, isend);
	const struct delta *delta;
	size_t j;
	const double *xj;
	const double *x;
	const size_t *ind;
	size_t dim = design2_tvar_dim(d);
	size_t nz;
	ptrdiff_t iz;

	if (version_changed(v, &sm->version)) {
		sender_clear(sm, d);
	}

	double t0 = sm->tcur;
	double t = history_time(h);

	if (t0 == t)
		return;

	design2_get_tvar_matrix(d, isend, &x, &ind, &nz);

	for (delta = deltaset_head(ds); delta != NULL; delta = delta->prev) {
		if (delta->t < t0 || delta->t == -INFINITY)
			break;

		j = delta->i;
		iz = find_index(j, ind, nz);
		assert(iz >= 0);
		xj = x + iz * dim;

		//mlogitaug_set_x(&sm->mlogitaug, j, 0, dim, xj);
		mlogitaug_set_x(&sm->mlogitaug, j, xj);
	}

	sm->tcur = t;
}


static void sender_init(struct recv_model_sender *sm, size_t isend,
			const struct design2 *d, const struct recv_params *p,
			struct recv_model_cohort *cm,
			struct mlogitaug_work *work)
{
	assert(isend < design2_count1(d));
	struct history *h = design2_history(d);
	size_t dim = design2_tvar_dim(d);

	mlogitaug_init(&sm->mlogitaug, &cm->mlogit, dim, work);
	version_watch_init(&sm->version, history_version(h));
	sender_clear(sm, d);
	sender_set(sm, isend, d, p);
}


static void sender_deinit(struct recv_model_sender *sm, struct history *h)
{
	version_watch_deinit(&sm->version, history_version(h));
	mlogitaug_deinit(&sm->mlogitaug);
}


void recv_model_init(struct recv_model *m,
		     const struct recv_params *p,
		     struct design *r, struct design2 *d)
{
	assert(m);
	assert(r);
	assert(design_count(r));
	assert(d);
	assert(design2_count2(d) == design_count(r));

	size_t nrecv = design_count(r);
	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_tvar_dim(r);
	size_t dimd0 = design_trait_dim(r);
	size_t dimd1 = design_tvar_dim(r);
	size_t base_dim = dimr0 + dimr1 + dimd0;
	size_t aug_dim = dimd1;
	size_t dim = base_dim + aug_dim;
	size_t isend, nsend = design2_count1(d);
	size_t ic, nc = design2_cohort_count(d);

	m->recv = r;
	m->dyad = d;
	m->exclude_loops = 0;

	mlogit_work_init(&m->work, nrecv, base_dim);
	mlogitaug_work_init(&m->augwork, base_dim, aug_dim);

	m->params.recv.traits = xcalloc(dim, sizeof(double));
	m->params.recv.tvars = m->params.recv.traits + dimr0;
	m->params.dyad.traits = m->params.recv.tvars + dimr1;
	m->params.dyad.tvars = m->params.dyad.traits + dimd0;

	if (p && p->recv.traits) {
		memcpy(m->params.recv.traits, p->recv.traits, dimr0 * sizeof(double));
	}
	if (p && p->recv.tvars) {
		memcpy(m->params.recv.tvars, p->recv.tvars, dimr1 * sizeof(double));
	}
	if (p && p->dyad.traits) {
		memcpy(m->params.dyad.traits, p->dyad.traits, dimd0 * sizeof(double));
	}
	if (p && p->dyad.tvars) {
		memcpy(m->params.dyad.tvars, p->dyad.tvars, dimd1 * sizeof(double));
	}

	struct recv_model_cohort *cms = xcalloc(nc, sizeof(*cms));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cms[ic], ic, r, d, &m->params, &m->work);
	}
	m->cohort_models = cms;

	struct recv_model_sender *sms = xcalloc(nsend, sizeof(*sms));
	for (isend = 0; isend < nsend; isend++) {
		size_t c = recv_model_cohort(m, isend);
		struct recv_model_cohort *cm = &cms[c];
		sender_init(&sms[isend], isend, d, &m->params, cm, &m->augwork);
	}
	m->sender_models = sms;

	m->moments = 2;
}


void recv_model_deinit(struct recv_model *m)
{
	assert(m);

	const struct design *r = recv_model_design(m);
	struct history *h = design_history(r);

	struct recv_model_sender *sms = m->sender_models;
	size_t isend, nsend = recv_model_send_count(m);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&sms[isend], h);
	}
	free(sms);

	struct recv_model_cohort *cms = m->cohort_models;
	size_t ic, nc = cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cms[ic], h);
	}
	free(cms);

	free(m->params.recv.traits);
	mlogitaug_work_deinit(&m->augwork);
	mlogit_work_deinit(&m->work);
}


struct design *recv_model_design(const struct recv_model *m)
{
	return m->recv;
}


struct design2 *recv_model_design2(const struct recv_model *m)
{
	return m->dyad;
}


const struct recv_params *recv_model_params(const struct recv_model *m)
{
	assert(m);
	return &((struct recv_model *)m)->params;
}


size_t recv_model_send_count(const struct recv_model *m)
{
	const struct design2 *d = recv_model_design2(m);
	return design2_count1(d);
}

size_t recv_model_count(const struct recv_model *m)
{
	const struct design *r = recv_model_design(m);
	return design_count(r);
}

size_t recv_model_dim(const struct recv_model *m)
{
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	
	size_t dim = (design_trait_dim(r) + design_tvar_dim(r)
		      + design2_trait_dim(d) + design2_tvar_dim(d));
	return dim;
}


struct catdist1 *recv_model_dist(const struct recv_model *m, size_t isend)
{
	assert(isend < recv_model_send_count(m));
	struct mlogitaug *mlogitaug = recv_model_mlogit(m, isend);
	struct catdist1 *dist = mlogitaug_dist(mlogitaug);
	return dist;
}


struct mlogitaug *recv_model_mlogit(const struct recv_model *m, size_t i)
{
	assert(i < recv_model_send_count(m));

	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	size_t c = recv_model_cohort(m, i);

	cohort_update(&m->cohort_models[c], c, r, d);
	sender_update(&m->sender_models[i], i, d);

	struct recv_model_sender *sm = &m->sender_models[i];
	struct mlogitaug *mlogitaug = &sm->mlogitaug;
	return mlogitaug;
}


void recv_model_set_params(struct recv_model *m, const struct recv_params *p)
{
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_tvar_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dimd1 = design2_tvar_dim(d);
	size_t dim = dimr0 + dimr1 + dimd0 + dimd1;

	memset(m->params.recv.traits, 0, dim * sizeof(double));

	if (p && p->recv.traits) {
		memcpy(m->params.recv.traits, p->recv.traits, dimr0 * sizeof(double));
	}
	if (p && p->recv.tvars) {
		memcpy(m->params.recv.tvars, p->recv.tvars, dimr1 * sizeof(double));
	}
	if (p && p->dyad.traits) {
		memcpy(m->params.dyad.traits, p->dyad.traits, dimd0 * sizeof(double));
	}
	if (p && p->dyad.tvars) {
		memcpy(m->params.dyad.tvars, p->dyad.tvars, dimd1 * sizeof(double));
	}

	struct recv_model_cohort *cms = m->cohort_models;
	size_t ic, nc = cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		cohort_set(&cms[ic], ic, r, d, &m->params);
	}

	struct recv_model_sender *sms = m->sender_models;
	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sender_set(&sms[is], is, d, &m->params);
	}
}


void recv_model_set_exclude_loops(struct recv_model *m, int exclude_loops)
{
	if (m->exclude_loops == exclude_loops)
		return;

	struct recv_model_sender *sms = m->sender_models;
	size_t is, ns = recv_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sender_set_exclude_loops(&sms[is], is, exclude_loops);
	}

	m->exclude_loops = exclude_loops;
}


int recv_model_moments(const struct recv_model *m)
{
	return m->moments;
}


void recv_model_set_moments(struct recv_model *m, int k)
{
	if (m->moments == k)
		return;

	struct recv_model_cohort *cms = m->cohort_models;
	size_t ic, nc = cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		cohort_set_moments(&cms[ic], k);
	}

	m->moments = k;
}


void recv_params_init(struct recv_params *p, const struct design *r, const struct design2 *d)
{
	coefs_init(&p->recv, r);
	coefs2_init(&p->dyad, d);
}


void recv_params_deinit(struct recv_params *p)
{
	coefs2_deinit(&p->dyad);
	coefs_deinit(&p->recv);
}



