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


static void cohort_set(struct recv_model_cohort *cm, size_t c,
		       const struct frame *f,
		       const struct recv_coefs *coefs)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t isend = design2_cohort_rep(d, c);
	size_t nrecv = design_count(r);
	size_t dimr = design_dim(r);
	size_t dimr0 = design_trait_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dim = dimr + dimd0;
	const double *xr = design_traits(r);
	const double *xd = design2_traits(d, isend);

	double *x = xcalloc(nrecv * dim, sizeof(*x));

	if (dimr0 > 0) {
		assert(dim > 0);
		lapack_dlacpy(LA_COPY_ALL, dimr0, nrecv, xr, dimr0, x, dim);
	}

	if (dimd0 > 0) {
		assert(dim > 0);
		lapack_dlacpy(LA_COPY_ALL, dimd0, nrecv, xd, dimd0, x + dimr, dim);
	}

	mlogit_set_all_offset(&cm->mlogit, NULL);
	mlogit_set_coefs(&cm->mlogit, coefs->all);
	mlogit_set_all_x(&cm->mlogit, x);
	free(x);
	//mlogit_check(&cm->mlogit);
}


static void cohort_clear(struct recv_model_cohort *cm)
{
	mlogit_set_all_x(&cm->mlogit, NULL);
}


static void cohort_init(struct recv_model_cohort *cm,
			size_t c,
			const struct frame *f,
			const struct recv_coefs *coefs)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t nrecv = frame_recv_count(f);
	size_t dimr = design_dim(r);
	size_t dimd0 = design2_trait_dim(d);
	size_t dim = dimr + dimd0;

	mlogit_init(&cm->mlogit, nrecv, dim);
	cohort_set(cm, c, f, coefs);
}


static void cohort_deinit(struct recv_model_cohort *cm)
{
	mlogit_deinit(&cm->mlogit);
}


size_t recv_model_cohort(const struct recv_model *m, size_t isend)
{
	const struct frame *f = recv_model_frame(m);
	const struct design2 *d = frame_dyad_design(f);
	return design2_cohort(d, isend);
}


static void sender_clear(struct recv_model_sender *sm)
{
	mlogitaug_set_all_x(&sm->mlogitaug, NULL, NULL, 0);
}


static void sender_set(struct recv_model_sender *sm,
		       size_t isend,
		       const struct frame *f, const struct recv_coefs *coefs)
{
	const struct design2 *d = frame_dyad_design(f);
	const int has_loops = frame_has_loops(f);

	const double *x;
	const size_t *indx;
	size_t nzx;
	design2_tvars_get(d, isend, &x, &indx, &nzx);
	mlogitaug_set_all_x(&sm->mlogitaug, indx, x, nzx);

	if (!has_loops) {
		const double neginf = -INFINITY;
		mlogitaug_set_all_offset(&sm->mlogitaug, &isend, &neginf, 1);
	}

	mlogitaug_set_coefs(&sm->mlogitaug, coefs->dyad.tvars);
	//mlogitaug_check(&sm->mlogitaug);
}


static void sender_init(struct recv_model_sender *sm, size_t isend,
			const struct frame *f, const struct recv_coefs *coefs,
			const struct recv_model_cohort *cm)
{
	assert(isend < frame_send_count(f));

	const struct design2 *d = frame_dyad_design(f);
	size_t dim = design2_tvar_dim(d);

	mlogitaug_init(&sm->mlogitaug, &cm->mlogit, dim);
	sender_set(sm, isend, f, coefs);
}


static void sender_deinit(struct recv_model_sender *sm)
{
	mlogitaug_deinit(&sm->mlogitaug);
}


static void recv_model_dyad_update(void *udata, struct design2 *d, size_t isend,
				   size_t jrecv, const double *delta,
				   const size_t *ind, size_t nz)
{
	struct recv_model *m = udata;
	assert(isend < recv_model_send_count(m));
	assert(jrecv < recv_model_send_count(m));

	struct recv_model_sender *sm = &m->sender_models[isend];
	mlogitaug_inc_x(&sm->mlogitaug, jrecv, ind, delta, nz);
	//mlogitaug_check(&sm->mlogitaug);
	(void)d;
}


static void recv_model_recv_update(void *udata, struct design *r,
				   size_t jrecv, const double *delta,
				   const size_t *ind, size_t nz)
{
	struct recv_model *m = udata;
	assert(jrecv < recv_model_send_count(m));
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

	struct recv_model_cohort *cms = m->cohort_models;

	size_t ic, nc = recv_model_cohort_count(m);

	for (ic = 0; ic < nc; ic++) {
		struct recv_model_cohort *cm = &cms[ic];
		mlogit_inc_x(&cm->mlogit, jrecv, m->ind_buf, delta, nz);
	}
	//mlogitaug_check(&sm->mlogitaug);
	(void)r;
}


static void recv_model_dyad_clear(void *udata, struct design2 *d)
{
	struct recv_model *m = udata;
	struct recv_model_sender *sms = m->sender_models;
	size_t i, n = recv_model_send_count(m);

	for (i = 0; i < n; i++) {
		sender_clear(&sms[i]);
	}

	(void)d;
}


static void recv_model_recv_clear(void *udata, struct design *d)
{
	struct recv_model *m = udata;
	struct recv_model_cohort *cms = m->cohort_models;
	size_t ic, nc = recv_model_cohort_count(m);

	for (ic = 0; ic < nc; ic++) {
		cohort_clear(&cms[ic]);
	}

	(void)d;
}


static const struct design_callbacks RECV_CALLBACKS = {
	recv_model_recv_update,
	NULL,
	recv_model_recv_clear
};

static const struct design2_callbacks DYAD_CALLBACKS = {
	recv_model_dyad_update,
	NULL,
	recv_model_dyad_clear
};


void recv_model_init(struct recv_model *model,
		     struct frame *f,
		     const struct recv_coefs *coefs)
{
	assert(model);
	assert(f);
	assert(frame_recv_count(f) > 0);
	assert(!frame_has_loops(f) || frame_recv_count(f) > 1);

	struct design *r = frame_recv_design(f);
	struct design2 *d = frame_dyad_design(f);
	size_t isend, nsend = frame_send_count(f);
	size_t ic, nc = design2_cohort_count(d);

	model->frame = f;

	recv_coefs_init(&model->coefs, f);
	size_t coefs_len = model->coefs.dim * sizeof(*model->coefs.all);
	if (coefs) {
		memcpy(model->coefs.all, coefs->all, coefs_len);
	} else {
		memset(model->coefs.all, 0, coefs_len);
	}

	struct recv_model_cohort *cms = xcalloc(nc, sizeof(*cms));
	for (ic = 0; ic < nc; ic++) {
		cohort_init(&cms[ic], ic, f, &model->coefs);
	}
	model->cohort_models = cms;

	struct recv_model_sender *sms = xcalloc(nsend, sizeof(*sms));
	for (isend = 0; isend < nsend; isend++) {
		size_t c = recv_model_cohort(model, isend);
		const struct recv_model_cohort *cm = &cms[c];
		sender_init(&sms[isend], isend, f, &model->coefs, cm);
	}
	model->sender_models = sms;

	model->ind_buf = xmalloc(coefs->dim * sizeof(*model->ind_buf));

	design_add_observer(r, model, &RECV_CALLBACKS);
	design2_add_observer(d, model, &DYAD_CALLBACKS);
}


void recv_model_deinit(struct recv_model *model)
{
	assert(model);

	struct frame *f = model->frame;
	struct design *r = frame_recv_design(f);
	struct design2 *d = frame_dyad_design(f);

	design2_remove_observer(d, model);
	design_remove_observer(r, model);

	free(model->ind_buf);

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
	const struct design2 *d = frame_dyad_design(f);
	return design2_cohort_count(d);
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

struct catdist1 *recv_model_dist(const struct recv_model *m, size_t isend)
{
	assert(isend < recv_model_send_count(m));
	const struct recv_model_sender *sm = &m->sender_models[isend];
	struct catdist1 *dist = mlogitaug_dist(&sm->mlogitaug);
	catdist1_update_cache(dist);

	return dist;
}
