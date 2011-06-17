#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sloglik.h"
#include "util.h"

void recv_sloglik_init(struct recv_sloglik *ll, const struct model *model, ssize_t isend)
{
	assert(ll);
	assert(model);
	assert(0 <= isend && isend < design_send_count(model_design(model)));

	ssize_t n = model_receiver_count(model);
	ssize_t p = model_dim(model);
	const struct design *design = model_design(model);
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	
	ll->model = (struct model *)model;
	
	ll->f = 0.0;
	vector_init(&ll->grad, p);
	ll->grad_cached = false;
	
	ll->nsend = 0;
	svector_init(&ll->nrecv, n);
	vector_init(&ll->dxobs, dyn_dim);
	
	ll->gamma = 0.0;
	svector_init(&ll->dp, n);
	vector_init(&ll->dxbar, dyn_dim);
}

void recv_sloglik_deinit(struct recv_sloglik *ll)
{
	assert(ll);

	vector_deinit(&ll->dxbar);
	svector_deinit(&ll->dp);
	vector_deinit(&ll->dxobs);
	svector_deinit(&ll->nrecv);
	vector_deinit(&ll->grad);
}

void recv_sloglik_add(struct recv_sloglik *ll,
		     const struct frame *f, ssize_t *jrecv, ssize_t n)
{
	ssize_t isend = ll->isend;
	const struct recv_model *model = model_recv_model(ll->model, f, isend);
	ssize_t nreceiver = recv_model_count(model);
	struct svector *wt = svector_alloc(nreceiver);
	double ntot = ll->nsend + n;
	double scale1 = n / ntot;
	double scale0 = 1 - scale1;
	double lpbar = 0.0;
	ssize_t i;

	ll->grad_cached = false;

	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = recv_model_logprob(model, jrecv[i]);
		lpbar += (lp - lpbar) / (i + 1);

		//*svector_item_ptr(wt, jrecv[i]) += 1.0;

		// update number of receives
		*svector_item_ptr(&ll->nrecv, jrecv[i]) += 1.0;
	}

	// update log likelihood
	ll->f += scale1 * ((-lpbar) - ll->f);

	// update observed variable diffs
	frame_recv_dmuls(scale1 / n, TRANS_TRANS, f, model->isend,
			 wt, scale0, &ll->dxobs);
	ll->gamma += scale1 * (model->gamma - ll->gamma);

	svector_scale(&ll->dp, scale0);
	
	ssize_t nactive = array_count(&model->active);
	for (i = 0; i < nactive; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(&model->active, i);
		double dp = vector_item(&model->dp, i);
		double *dst = svector_item_ptr(&ll->dp, jrecv);
		*dst += scale1 * dp;
	}
	//svector_axpys(scale1, &model->dp, &ll->dp);

	vector_scale(&ll->dxbar, scale0);
	vector_axpy(scale1, &model->mean_dx, &ll->dxbar);

	// update number of sends
	ll->nsend += n;

	svector_free(wt);
}

double recv_sloglik_value(const struct recv_sloglik *ll)
{
	assert(ll);
	return (-ll->nsend) * ll->f;
}

/*
 *         g = [ sum{gamma[t,i]} * xbar[0,i]
 *               + ( X[0,i])^T * sum{dP[t,i]}
 *               + sum{dxbar[t,i]} ]
 *             -
 *             [ (X[0,i])^T n[i] + sum{dx[t,i,j]} ]
 */
static void
recv_sloglik_axpy_grad_nocache(double alpha, const struct recv_sloglik *ll, struct vector *y)
{
	double scale = (-ll->nsend) * alpha;
	const struct recv_model *rm = model_recv_model(ll->model, NULL, ll->isend);
	const struct vector *xbar0 = recv_model_mean0(rm);

	const struct design *design = model_design(ll->model);
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);
	struct vector ysub = vector_slice(y, off, dim);


	// sum{gamma[t,i]} * xbar[0,i]
	vector_axpy(scale * ll->gamma, xbar0, y);

	// (X[0,i])^T * sum{dP[t,i]}
	design_recv_muls0(scale, TRANS_TRANS,
			  ll->model->design, ll->isend, &ll->dp, 1.0,
			  y);

	// sum{dxbar[t,i]}
	vector_axpy(scale, &ll->dxbar, &ysub);

	// - (X[0,i])^T n[i]
	design_recv_muls0(-scale / ll->nsend, TRANS_TRANS,
			  ll->model->design, ll->isend, &ll->nrecv,
			  1.0, y);

	// -sum{dx[t,i,j]}
	vector_axpy(-scale, &ll->dxobs, &ysub);
}

static void
recv_sloglik_cache_grad(struct recv_sloglik *ll)
{
	vector_fill(&ll->grad, 0.0);
	recv_sloglik_axpy_grad_nocache(1.0, ll, &ll->grad);
	ll->grad_cached = true;
}

void recv_sloglik_axpy_grad(double alpha, const struct recv_sloglik *ll, struct vector *y)
{
	assert(ll);
	assert(y);

	if (!ll->grad_cached) {
		recv_sloglik_cache_grad((struct recv_sloglik *)ll);
	}
	
	vector_axpy(alpha, &ll->grad, y);
}
