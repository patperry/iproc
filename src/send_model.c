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
#include "send_model.h"



size_t send_params_dim(const struct design *d)
{
	return design_dim(d);
}


void send_params_init(struct send_params *p, const struct design *d)
{
	size_t dim = send_params_dim(d);
	double *data = xmalloc(dim * sizeof(double));
	send_params_init_view(p, d, data);
	p->owner = 1;
}


void send_params_init_view(struct send_params *p, const struct design *d,
			  const double *data)
{
	p->all = (double *)data;
	p->dim = 0;

	coefs_init_view(&p->coefs, d, p->all + p->dim);
	p->dim += p->coefs.dim;

	p->owner = 0;
}


void send_params_deinit(struct send_params *p)
{
	if (p->owner)
		free(p->all);
}


static void base_set(struct mlogit *base,
		     const struct design *d,
		     const struct send_params *p)
{
	const double *x = design_trait_matrix(d);
	const double *beta = p->coefs.all;

	mlogit_set_all_offset(base, NULL);
	mlogit_set_all_x(base, x);
	mlogit_set_coefs(base, beta);
	//mlogit_check(&cm->mlogit);
}


static void cohort_init(struct send_model_cohort *cm,
			size_t c,
			const struct frame *f,
			const struct send_params *coefs,
			struct mlogit_work *work)
{
	const struct design *r = frame_send_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t nsend = frame_send_count(f);
	size_t dimr = design_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dim = dimr + dimd0;

	mlogit_init(&cm->mlogit, nsend, dim, work);
	cohort_set(cm, c, f, coefs);
}


static void cohort_deinit(struct send_model_cohort *cm)
{
	mlogit_deinit(&cm->mlogit);
}


size_t send_model_cohort(const struct send_model *m, size_t isend)
{
	const struct frame *f = send_model_frame(m);
	const struct design2 *d = frame_dyad_design(f);
	return design2_cohort(d, isend);
}


static void sender_clear(struct send_model_sender *sm)
{
	mlogitaug_set_all_x(&sm->mlogitaug, NULL, NULL, 0);
}


static void sender_set(struct send_model_sender *sm,
		       size_t isend,
		       const struct frame *f, const struct send_params *coefs)
{
	const struct design2 *d = frame_dyad_design(f);
	const int has_loops = frame_has_loops(f);

	const double *x;
	const size_t *indx;
	size_t nzx;
	design2_tvars_get_all(d, isend, &x, &indx, &nzx);
	mlogitaug_set_all_x(&sm->mlogitaug, indx, x, nzx);

	if (!has_loops) {
		const double neginf = -INFINITY;
		mlogitaug_set_all_offset(&sm->mlogitaug, &isend, &neginf, 1);
	}

	mlogitaug_set_coefs(&sm->mlogitaug, coefs->dyad.tvars);
	//mlogitaug_check(&sm->mlogitaug);
}


static void sender_init(struct send_model_sender *sm, size_t isend,
			const struct frame *f, const struct send_params *coefs,
			struct send_model_cohort *cm,
			struct mlogitaug_work *work)
{
	assert(isend < frame_send_count(f));

	const struct design2 *d = frame_dyad_design(f);
	size_t dim = design2_tvar_dim(d);

	mlogitaug_init(&sm->mlogitaug, &cm->mlogit, dim, work);
	sender_set(sm, isend, f, coefs);
}


static void sender_deinit(struct send_model_sender *sm)
{
	mlogitaug_deinit(&sm->mlogitaug);
}


static void send_model_dyad_update(void *udata, struct design2 *d, size_t isend,
				   size_t jsend, const double *delta,
				   const size_t *ind, size_t nz)
{
	struct send_model *m = udata;
	assert(isend < send_model_send_count(m));
	assert(jsend < send_model_send_count(m));

	struct send_model_sender *sm = &m->sender_models[isend];
	mlogitaug_inc_x(&sm->mlogitaug, jsend, ind, delta, nz);
	//mlogitaug_check(&sm->mlogitaug);
	(void)d;
}


static void send_model_send_update(void *udata, struct design *r,
				   size_t jsend, const double *delta,
				   const size_t *ind, size_t nz)
{
	struct send_model *m = udata;
	assert(jsend < send_model_send_count(m));
	assert(nz <= design_tvar_dim(r));

	size_t off = design_trait_dim(r);
	size_t iz;

	if (ind) {
		for (iz = 0; iz < nz; iz++) {
			m->ind_buf[iz] = off + ind[iz];
		}
	} else {
		assert(nz == design_tvar_dim(r));
		for (iz = 0; iz < nz; iz++) {
			m->ind_buf[iz] = off + iz;
		}
	}

	struct send_model_cohort *cms = m->cohort_models;

	size_t ic, nc = send_model_cohort_count(m);

	for (ic = 0; ic < nc; ic++) {
		struct send_model_cohort *cm = &cms[ic];
		mlogit_inc_x(&cm->mlogit, jsend, m->ind_buf, delta, nz);
	}
	//mlogitaug_check(&sm->mlogitaug);
	(void)r;
}


static void send_model_dyad_clear(void *udata, struct design2 *d)
{
	struct send_model *m = udata;
	struct send_model_sender *sms = m->sender_models;
	size_t i, n = send_model_send_count(m);

	for (i = 0; i < n; i++) {
		sender_clear(&sms[i]);
	}

	(void)d;
}


static void send_model_send_clear(void *udata, struct design *d)
{
	struct send_model *m = udata;
	struct frame *f = m->frame;
	struct send_model_cohort *cms = m->cohort_models;
	size_t ic, nc = send_model_cohort_count(m);

	for (ic = 0; ic < nc; ic++) {
		cohort_clear(&cms[ic], ic, f);
	}

	(void)d;
}


static const struct design_callbacks RECV_CALLBACKS = {
	send_model_send_update,
	NULL,
	send_model_send_clear
};

static const struct design2_callbacks DYAD_CALLBACKS = {
	send_model_dyad_update,
	NULL,
	send_model_dyad_clear
};


void send_model_init(struct send_model *model,
		     struct frame *f,
		     const struct send_params *coefs)
{
	assert(model);
	assert(f);
	assert(frame_send_count(f) > 0);
	assert(!frame_has_loops(f) || frame_send_count(f) > 1);

	struct design *r = frame_send_design(f);
	struct design2 *d = frame_dyad_design(f);
	size_t nsend = frame_send_count(f);
	size_t base_dim = design_dim(r) + design2_trait_dim(d);
	size_t aug_dim = design2_tvar_dim(d);
	size_t isend, nsend = frame_send_count(f);
	size_t ic, nc = design2_cohort_count(d);

	model->frame = f;

	mlogit_work_init(&model->work, nsend, base_dim);
	mlogitaug_work_init(&model->augwork, base_dim, aug_dim);

	send_params_init(&model->coefs, f);
	size_t coefs_len = model->coefs.dim * sizeof(*model->coefs.all);
	if (coefs) {
		memcpy(model->coefs.all, coefs->all, coefs_len);
	} else {
		memset(model->coefs.all, 0, coefs_len);
	}

	struct send_model_cohort *cms = xcalloc(nc, sizeof(*cms));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cms[ic], ic, f, &model->coefs, &model->work);
	}
	model->cohort_models = cms;

	struct send_model_sender *sms = xcalloc(nsend, sizeof(*sms));
	for (isend = 0; isend < nsend; isend++) {
		size_t c = send_model_cohort(model, isend);
		struct send_model_cohort *cm = &cms[c];
		sender_init(&sms[isend], isend, f, &model->coefs, cm, &model->augwork);
	}
	model->sender_models = sms;

	model->ind_buf = xmalloc(model->coefs.dim * sizeof(*model->ind_buf));

	model->moments = 2;

	design_add_observer(r, model, &RECV_CALLBACKS);
	design2_add_observer(d, model, &DYAD_CALLBACKS);
}


void send_model_deinit(struct send_model *model)
{
	assert(model);

	struct frame *f = model->frame;
	struct design *r = frame_send_design(f);
	struct design2 *d = frame_dyad_design(f);

	design2_remove_observer(d, model);
	design_remove_observer(r, model);

	free(model->ind_buf);

	struct send_model_sender *sms = model->sender_models;
	size_t isend, nsend = send_model_send_count(model);
	for (isend = 0; isend < nsend; isend++) {
		sender_deinit(&sms[isend]);
	}
	free(sms);

	struct send_model_cohort *cms = model->cohort_models;
	size_t ic, nc = send_model_cohort_count(model);
	for (ic = 0; ic < nc; ic++) {
		cohort_deinit(&cms[ic]);
	}
	free(cms);

	send_params_deinit(&model->coefs);
	mlogitaug_work_deinit(&model->augwork);
	mlogit_work_deinit(&model->work);
}


struct frame *send_model_frame(const struct send_model *model)
{
	assert(model);
	return model->frame;
}

const struct send_params *send_model_coefs(const struct send_model *model)
{
	assert(model);
	return &((struct send_model *)model)->coefs;
}

size_t send_model_send_count(const struct send_model *model)
{
	assert(model);
	return frame_send_count(model->frame);
}

size_t send_model_cohort_count(const struct send_model *model)
{
	assert(model);
	const struct frame *f = send_model_frame(model);
	const struct design2 *d = frame_dyad_design(f);
	return design2_cohort_count(d);
}

size_t send_model_count(const struct send_model *model)
{
	assert(model);
	const struct frame *f = send_model_frame(model);
	const struct design *r = frame_send_design(f);
	return design_count(r);
}

size_t send_model_dim(const struct send_model *m)
{
	assert(m);
	const struct frame *f = send_model_frame(m);
	const struct design *r = frame_send_design(f);
	const struct design2 *d = frame_dyad_design(f);

	size_t dim = (design_trait_dim(r) + design_tvar_dim(r)
		      + design2_trait_dim(d) + design2_tvar_dim(d));
	return dim;
}

struct catdist1 *send_model_dist(const struct send_model *m, size_t isend)
{
	assert(isend < send_model_send_count(m));
	struct mlogitaug *mlogitaug = send_model_mlogit(m, isend);
	struct catdist1 *dist = mlogitaug_dist(mlogitaug);
	return dist;
}

struct mlogitaug *send_model_mlogit(const struct send_model *m, size_t isend)
{
	assert(isend < send_model_send_count(m));
	struct send_model_sender *sm = &m->sender_models[isend];
	struct mlogitaug *mlogitaug = &sm->mlogitaug;
	return mlogitaug;
}

void send_model_set_params(struct send_model *m, const struct send_params *coefs)
{
	const struct frame *f = send_model_frame(m);
	size_t dim = send_model_dim(m);

	if (coefs) {
		memcpy(m->coefs.all, coefs->all, dim * sizeof(*m->coefs.all));
	} else {
		memset(m->coefs.all, 0, dim * sizeof(*m->coefs.all));
	}

	struct send_model_cohort *cms = m->cohort_models;
	size_t ic, nc = send_model_cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		cohort_set(&cms[ic], ic, f, &m->coefs);
	}

	struct send_model_sender *sms = m->sender_models;
	size_t is, ns = send_model_send_count(m);
	for (is = 0; is < ns; is++) {
		sender_set(&sms[is], is, f, &m->coefs);
	}
}

int send_model_moments(const struct send_model *m)
{
	return m->moments;
}

void send_model_set_moments(struct send_model *m, int k)
{
	struct send_model_cohort *cms = m->cohort_models;
	size_t ic, nc = send_model_cohort_count(m);
	for (ic = 0; ic < nc; ic++) {
		cohort_set_moments(&cms[ic], k);
	}

	m->moments = k;
}
