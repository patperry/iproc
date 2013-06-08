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


void send_model_update(struct send_model *m)
{
	const struct design *d = send_model_design(m);
	const struct history *h = design_history(d);
	const struct version *v = history_version(h);
	const struct deltaset *ds = design_changes(d);
	const struct delta *delta;
	double t0 = m->tcur;
	double t = history_time(h);
	size_t i;
	const double *xi;
	const double *x;
	const size_t *ind;
	size_t dim = design_tvar_dim(d);
	size_t nz;
	ptrdiff_t iz;

	assert(t0 <= t);
	if (t0 == t && !version_changed(v, &m->history_version))
		return;

	design_get_tvar_matrix(d, &x, &ind, &nz);
	version_watch_set(&m->history_version, v);

	for (delta = deltaset_head(ds); delta != NULL; delta = delta->prev) {
		if (delta->t < t0 || delta->t == -INFINITY)
			break;

		i = delta->i;
		iz = find_index(i, ind, nz);
		assert(iz >= 0);
		xi = x + iz * dim;

		mlogitaug_set_x(&m->aug, i, xi);
	}

	m->tcur = t;
}

//static void sender_clear(struct send_model_sender *sm)
//{
//	mlogitaug_set_all_x(&sm->mlogitaug, NULL, NULL, 0);
//}



void send_model_init(struct send_model *m, const struct send_params *p,
		     struct design *d)
{
	assert(design_count(d));

	const struct history *h = design_history(d);
	struct version *v = history_version(h);
	size_t count = design_count(d);
	size_t base_dim = design_trait_dim(d);
	const double *trait_x = design_trait_matrix(d);
	size_t aug_dim = design_tvar_dim(d);
	size_t dim = base_dim + aug_dim;

	m->design = d;
	version_watch_init(&m->history_version, v);

	m->params.coefs.traits = xcalloc(dim, sizeof(double));
	m->params.coefs.tvars = m->params.coefs.traits + base_dim;

	if (p) {
		memcpy(m->params.coefs.traits, p->coefs.traits,
		       base_dim * sizeof(double));
		memcpy(m->params.coefs.tvars, p->coefs.tvars,
		       aug_dim * sizeof(double));
	}

	mlogit_work_init(&m->base_work, count, base_dim);
	mlogit_init(&m->base, count, base_dim, &m->base_work);
	mlogit_set_all_x(&m->base, trait_x);
	mlogit_set_coefs(&m->base, m->params.coefs.traits);

	mlogitaug_work_init(&m->aug_work, base_dim, aug_dim);
	mlogitaug_init(&m->aug, &m->base, aug_dim, &m->aug_work);
	mlogitaug_set_coefs(&m->aug, m->params.coefs.tvars);

	m->tcur = -INFINITY;
}


void send_model_deinit(struct send_model *m)
{
	const struct design *d = send_model_design(m);
	const struct history *h = design_history(d);
	struct version *v = history_version(h);

	assert(m);

	mlogitaug_deinit(&m->aug);
	mlogitaug_work_deinit(&m->aug_work);
	mlogit_deinit(&m->base);
	mlogit_work_deinit(&m->base_work);
	free(m->params.coefs.traits);
	version_watch_deinit(&m->history_version, v);
}


struct design *send_model_design(const struct send_model *m)
{
	return m->design;
}

const struct send_params *send_model_params(const struct send_model *m)
{
	return &((struct send_model *)m)->params;
}

size_t send_model_count(const struct send_model *m)
{
	return design_count(m->design);
}

size_t send_model_dim(const struct send_model *m)
{
	return design_dim(m->design);
}

struct catdist1 *send_model_dist(const struct send_model *m)
{
	const struct mlogitaug *aug = send_model_mlogit(m);
	struct catdist1 *dist = mlogitaug_dist(aug);
	return dist;
}

struct mlogitaug *send_model_mlogit(const struct send_model *m)
{
	send_model_update((struct send_model *)m);
	return (struct mlogitaug *)&m->aug;
}

void send_model_set_params(struct send_model *m, const struct send_params *p)
{
	const struct design *d = send_model_design(m);
	size_t base_dim = design_trait_dim(d);
	size_t aug_dim = design_tvar_dim(d);

	if (p) {
		memcpy(m->params.coefs.traits, p->coefs.traits,
		       base_dim * sizeof(double));
		memcpy(m->params.coefs.tvars, p->coefs.tvars,
		       aug_dim * sizeof(double));
	} else {
		memset(m->params.coefs.traits, 0, base_dim * sizeof(double));
		memset(m->params.coefs.tvars, 0, aug_dim * sizeof(double));
	}

	mlogit_set_coefs(&m->base, m->params.coefs.traits);
	mlogitaug_set_coefs(&m->aug, m->params.coefs.tvars);
}


int send_model_moments(const struct send_model *m)
{
	return mlogit_moments(&m->base);
}

void send_model_set_moments(struct send_model *m, int k)
{
	mlogit_set_moments(&m->base, k);
}
