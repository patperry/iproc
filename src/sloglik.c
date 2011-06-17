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
	
	ll->model = (struct model *)model;
	
	ll->f = 0.0;
	vector_init(&ll->grad, p);
	ll->grad_cached = false;
	
	ll->nsend = 0;
	svector_init(&ll->nrecv, n);
	svector_init(&ll->dxobs, p);
	
	ll->gamma = 0.0;
	svector_init(&ll->dp, n);
	svector_init(&ll->dxbar, p);
}

void recv_sloglik_deinit(struct recv_sloglik *ll)
{
	assert(ll);

	svector_deinit(&ll->dxbar);
	svector_deinit(&ll->dp);
	svector_deinit(&ll->dxobs);
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

		*svector_item_ptr(wt, jrecv[i]) += 1.0;

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
	svector_axpys(scale1, &model->dp, &ll->dp);

	svector_scale(&ll->dxbar, scale0);
	svector_axpys(scale1, &model->dxbar, &ll->dxbar);

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

	// sum{gamma[t,i]} * xbar[0,i]
	vector_axpy(scale * ll->gamma, xbar0, y);

	// (X[0,i])^T * sum{dP[t,i]}
	design_recv_muls0(scale, TRANS_TRANS,
			  ll->model->design, ll->isend, &ll->dp, 1.0,
			  y);

	// sum{dxbar[t,i]}
	svector_axpy(scale, &ll->dxbar, y);

	// - (X[0,i])^T n[i]
	design_recv_muls0(-scale / ll->nsend, TRANS_TRANS,
			  ll->model->design, ll->isend, &ll->nrecv,
			  1.0, y);

	// -sum{dx[t,i,j]}
	svector_axpy(-scale, &ll->dxobs, y);
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
