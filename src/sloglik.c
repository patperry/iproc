#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sloglik.h"
#include "util.h"


static void mean_init(struct recv_sloglik_mean *mean,
		      const struct design *design)
{
	assert(mean);
	assert(design);
	
	ssize_t dyn_dim = design_recv_dyn_dim(design);
	
	mean->gamma = 0.0;
	vector_init(&mean->dp, 0);
	vector_init(&mean->mean_dx, dyn_dim);
}

static void mean_deinit(struct recv_sloglik_mean *mean)
{
	assert(mean);

	vector_deinit(&mean->mean_dx);
	vector_deinit(&mean->dp);
}

static void mean_clear(struct recv_sloglik_mean *mean)
{
	assert(mean);
	mean->gamma = 0.0;
	vector_reinit(&mean->dp, 0);
	vector_fill(&mean->mean_dx, 0.0);
}

static ssize_t mean_active_count(const struct recv_sloglik_mean *mean)
{
	assert(mean);
	return vector_dim(&mean->dp);
}

static void mean_insert_active(const struct recv_sloglik_mean *mean,
			       ssize_t index)
{
	assert(mean);
	assert(0 <= index && index <= mean_active_count(mean));
	
	ssize_t n0 = mean_active_count(mean);
	ssize_t n1 = n0 + 1;
	struct vector dp, src, dst;
	vector_init(&dp, n1);
	
	src = vector_slice(&mean->dp, 0, index);
	dst = vector_slice(&dp, 0, index);
	vector_assign_copy(&dst, &src);
	
	src = vector_slice(&mean->dp, index, n0 - index);
	dst = vector_slice(&dp, index + 1, n0 - index);
	vector_assign_copy(&dst, &src);	
}

static void mean_grow_active(const struct recv_sloglik_mean *mean,
			     const struct array *active0,
			     const struct array *active1)
{
	assert(mean);
	assert(active0);
	assert(active1);
	assert(mean_active_count(mean) == array_count(active0));
	assert(array_count(active0) <= array_count(active1));
	
	if (array_count(active0) == array_count(active1))
		return;
	
	const ssize_t n0 = array_count(active0);
	const ssize_t n1 = array_count(active1);
	const ssize_t *begin0 = array_to_ptr(active0);
	const ssize_t *end0 = begin0 + n0;
	const ssize_t *begin1 = array_to_ptr(active1);
	const ssize_t *end1 = begin1 + n1;
	const ssize_t *i0, *i1;
	
	for (i0 = begin0, i1 = begin1; i1 < end1; i1++) {
		if (i0 < end0 && *i0 == *i1) {
			i0++;
		} else {
			assert(i0 == end0 || *i1 < *i0);
			mean_insert_active(mean, i1 - begin1);
		}
	}
	assert(i0 == end0);
}

static void mean_set(struct recv_sloglik_mean *mean,
		     const struct recv_model *model,
		     const struct frame *f)
{
	assert(mean);
	assert(model);
	assert(f);
	
	const ssize_t isend = model->isend;
	const struct array *active = &model->active;
	const ssize_t n = array_count(active);
	ssize_t i;

	const double gamma = model->gamma;
	mean->gamma = gamma;
	
	if (vector_dim(&mean->dp) != n) {
		vector_reinit(&mean->dp, n);
	}
	vector_fill(&mean->mean_dx, 0.0);

	for (i = 0; i < n; i++) {
		ssize_t jrecv = *(ssize_t *)array_item(active, i);
		
		/* mean */		
		double p = recv_model_prob(model, jrecv);
		double p0 = recv_model_prob0(model, jrecv);
		double dp = p - gamma * p0;
		vector_set_item(&mean->dp, i, dp);
		
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		vector_axpy(p, dx, &mean->mean_dx);
	}
}

static void mean_scale(struct recv_sloglik_mean *mean, double scale)
{
	assert(mean);
	
	mean->gamma *= scale;
	vector_scale(&mean->dp, scale);
	vector_scale(&mean->mean_dx, scale);
}

static void mean_axpy(double alpha,
		      const struct recv_sloglik_mean *mean0,
		      struct recv_sloglik_mean *mean1)
{
	assert(mean0);
	assert(mean1);	
	assert(mean_active_count(mean0) == mean_active_count(mean1));

	mean1->gamma += alpha * (mean0->gamma);
	vector_axpy(alpha, &mean0->dp, &mean1->dp);
	vector_axpy(alpha, &mean0->mean_dx, &mean1->mean_dx);
}

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
	ll->isend = isend;
	
	ll->f = 0.0;
	vector_init(&ll->grad, p);
	ll->grad_cached = false;
	
	ll->nsend = 0;
	svector_init(&ll->nrecv, n);
	vector_init(&ll->dxobs, dyn_dim);
	
	ll->dev = 0.0;
	ll->dev_last = 0.0;
	array_init(&ll->active, sizeof(ssize_t));
	mean_init(&ll->mean_last, design);
	mean_init(&ll->mean, design);
	
	/* deprecated */
	ll->gamma = 0.0;
	svector_init(&ll->dp, n);
	vector_init(&ll->dxbar, dyn_dim);
}

void recv_sloglik_deinit(struct recv_sloglik *ll)
{
	assert(ll);

	mean_deinit(&ll->mean);
	mean_deinit(&ll->mean_last);
	array_deinit(&ll->active);
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
	double ntot = ll->nsend + n;
	double scale1 = n / ntot;
	//double scale0 = 1 - scale1;
	ssize_t i;

	ll->dev_last = 0.0;
	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);
		
		double lp = recv_model_logprob(model, jrecv[i]);
		ll->dev_last += -2 * lp;
	}
	
	ll->dev += scale1 * (ll->dev_last / n - ll->dev);

	ll->nsend = ntot;
	
	/*
	 struct svector *wt = svector_alloc(nreceiver);
	 double lpbar = 0.0;


	
	ll->grad_cached = false;

	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = recv_model_logprob(model, jrecv[i]);
		lpbar += (lp - lpbar) / (i + 1);

		// *svector_item_ptr(wt, jrecv[i]) += 1.0;

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

	svector_free(wt);*/
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

ssize_t recv_sloglik_count(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->nsend;
}

double recv_sloglik_mean_dev(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->dev;
}

double recv_sloglik_last_dev(const struct recv_sloglik *sll)
{
	assert(sll);
	return sll->dev_last;
}
