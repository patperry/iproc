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



//static void sender_clear(struct send_model_sender *sm)
//{
//	mlogitaug_set_all_x(&sm->mlogitaug, NULL, NULL, 0);
//}



void send_model_init(struct send_model *model,
		     struct design *d,
		     const struct send_params *p)
{
	assert(design_count(d));

	size_t count = design_count(d);
	size_t base_dim = design_trait_dim(d);
	const double *trait_x = design_trait_matrix(d);
	size_t aug_dim = design_tvar_dim(d);
	
	model->design = d;

	send_params_init(&model->params, d);
	size_t len = model->params.dim * sizeof(double);
	if (len) {
		memcpy(model->params.all, p->all, len);
	} else {
		memset(model->params.all, 0, len);
	}

	mlogit_work_init(&model->base_work, count, base_dim);
	mlogit_init(&model->base, count, base_dim, &model->base_work);
	mlogit_set_all_x(&model->base, trait_x);
	mlogit_set_coefs(&model->base, model->params.coefs.traits);

	mlogitaug_work_init(&model->aug_work, base_dim, aug_dim);
	mlogitaug_init(&model->aug, &model->base, aug_dim, &model->aug_work);
	mlogitaug_set_coefs(&model->aug, model->params.coefs.tvars);
}


void send_model_deinit(struct send_model *model)
{
	assert(model);

	mlogitaug_deinit(&model->aug);
	mlogitaug_work_deinit(&model->aug_work);
	mlogit_deinit(&model->base);
	mlogit_work_deinit(&model->base_work);
	send_params_deinit(&model->params);
}


struct design *send_model_design(const struct send_model *model)
{
	return model->design;
}

const struct send_params *send_model_params(const struct send_model *model)
{
	return &((struct send_model *)model)->params;
}

size_t send_model_count(const struct send_model *model)
{
	return design_count(model->design);
}

size_t send_model_dim(const struct send_model *m)
{
	return design_dim(m->design);
}

struct catdist1 *send_model_dist(const struct send_model *m)
{
	struct catdist1 *dist = mlogitaug_dist(&m->aug);
	return dist;
}

struct mlogitaug *send_model_mlogit(const struct send_model *m)
{
	return (struct mlogitaug *)&m->aug;
}

void send_model_set_params(struct send_model *m, const struct send_params *p)
{
	size_t dim = send_model_dim(m);

	if (p) {
		memcpy(m->params.all, p->all, dim * sizeof(double));
	} else {
		memset(m->params.all, 0, dim * sizeof(double));
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
