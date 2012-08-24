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
	const struct design *s = frame_send_design(f);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t isend = design_cohort_rep(s, c);

	design_traits_mul(1.0, r, coefs->recv.traits, 0.0, cm->eta0);
	design2_traits_mul(1.0, d, isend, coefs->dyad.traits, 1.0, cm->eta0);

	mlogit_set_all_offset(&cm->mlogit, cm->eta0);
	mlogit_set_coefs(&cm->mlogit, coefs->recv.tvars);
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
	size_t nrecv = frame_recv_count(f);
	size_t dim = coefs->recv.dim + coefs->dyad.dim;

	mlogit_init(&cm->mlogit, nrecv, dim);
	cm->eta0 = xmalloc(nrecv * sizeof(double));
	cohort_set(cm, c, f, coefs);
}


static void cohort_deinit(struct recv_model_cohort *cm)
{
	free(cm->eta0);
	mlogit_deinit(&cm->mlogit);
}


size_t recv_model_cohort(const struct recv_model *m, size_t isend)
{
	const struct frame *f = recv_model_frame(m);
	const struct design *s = frame_send_design(f);
	return design_cohort(s, isend);
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


static void model_clear(struct recv_model *m)
{
	struct recv_model_sender *sms = m->sender_models;
	size_t i, n = recv_model_send_count(m);

	for (i = 0; i < n; i++) {
		sender_clear(&sms[i]);
	}

	struct recv_model_cohort *cms = m->cohort_models;
	size_t c, nc = recv_model_cohort_count(m);

	for (c = 0; c < nc; c++) {
		cohort_clear(&cms[c]);
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
		size_t c = recv_model_cohort(model, isend);
		const struct recv_model_cohort *cm = &cms[c];
		sender_init(&sms[isend], isend, f, &model->coefs, cm);
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


struct catdist1 *recv_model_dist(const struct recv_model *m, size_t isend)
{
	assert(isend < recv_model_send_count(m));
	const struct recv_model_sender *sm = &m->sender_models[isend];
	struct catdist1 *dist = mlogitaug_dist(&sm->mlogitaug);
	return dist;
}


#if 0

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
		sender_set(&senders[isend], isend, f, coefs);
	}
}

struct recv_model_sender *recv_model_send(const struct recv_model *m,
					  size_t isend)
{
	return &m->sender_models[isend];
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

double recv_model_invgrow(const struct recv_model *m, size_t isend)
{
	assert(m);
	assert(isend < recv_model_send_count(m));

	const struct recv_model_sender *rm = recv_model_send(m, isend);
	return rm->gamma;
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

#endif
