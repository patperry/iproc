#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "sloglik.h"

static void iproc_sloglik_free(iproc_sloglik * sll)
{
	if (sll) {
		svector_free(sll->dxbar);
		svector_free(sll->dp);
		svector_free(sll->dxobs);
		svector_free(sll->nrecv);
		vector_free(sll->grad);
		iproc_model_unref(sll->model);
		free(sll);
	}
}

iproc_sloglik *iproc_sloglik_new(iproc_model * model, ssize_t isend)
{
	iproc_sloglik *sll = malloc(sizeof(*sll));
	if (!sll)
		return NULL;

	ssize_t n = iproc_model_nreceiver(model);
	ssize_t p = iproc_model_dim(model);

	sll->model = iproc_model_ref(model);
	sll->isend = isend;

	sll->f = 0.0;
	sll->grad = vector_alloc(p);
	sll->grad_cached = false;

	sll->nsend = 0;
	sll->nrecv = svector_alloc(n);
	sll->dxobs = svector_alloc(p);

	sll->gamma = 0.0;
	sll->dp = svector_alloc(n);
	sll->dxbar = svector_alloc(p);

	refcount_init(&sll->refcount);

	if (!(sll->grad && sll->nrecv && sll->dxobs && sll->dxbar)) {
		iproc_sloglik_free(sll);
		sll = NULL;
	}

	return sll;
}

iproc_sloglik *iproc_sloglik_ref(iproc_sloglik * sll)
{
	if (sll) {
		refcount_get(&sll->refcount);
	}
	return sll;
}

static void iproc_sloglik_release(struct refcount *refcount)
{
	iproc_sloglik *sll = container_of(refcount, iproc_sloglik, refcount);
	iproc_sloglik_free(sll);
}

void iproc_sloglik_unref(iproc_sloglik * sll)
{
	if (sll) {
		refcount_put(&sll->refcount, iproc_sloglik_release);
	}
}

void iproc_sloglik_insert(iproc_sloglik * sll,
			  const struct frame *f, ssize_t jrecv)
{
	iproc_sloglik_insertm(sll, f, &jrecv, 1);
}

void iproc_sloglik_insertm(iproc_sloglik * sll,
			   const struct frame *f, ssize_t *jrecv, ssize_t n)
{
	ssize_t nreceiver = iproc_model_nreceiver(sll->model);
	iproc_model_ctx *ctx = iproc_model_ctx_new(sll->model, f, sll->isend);
	struct svector *wt = svector_alloc(nreceiver);
	double ntot = sll->nsend + n;
	double scale1 = n / ntot;
	double scale0 = 1 - scale1;
	double lpbar = 0.0;
	ssize_t i;

	sll->grad_cached = false;

	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = iproc_model_ctx_logprob(ctx, jrecv[i]);
		lpbar += (lp - lpbar) / (i + 1);

		*svector_at(wt, jrecv[i]) += 1.0;

		// update number of receives
		*svector_at(sll->nrecv, jrecv[i]) += 1.0;
	}

	// update log likelihood
	sll->f += scale1 * ((-lpbar) - sll->f);

	// update observed variable diffs
	frame_dmuls(scale1 / n, TRANS_TRANS, f, sll->isend,
		    wt, scale0, sll->dxobs);
	sll->gamma += scale1 * (ctx->gamma - sll->gamma);

	svector_scale(sll->dp, scale0);
	svector_axpys(scale1, ctx->dp, sll->dp);

	svector_scale(sll->dxbar, scale0);
	svector_axpys(scale1, ctx->dxbar, sll->dxbar);

	// update number of sends
	sll->nsend += n;

	svector_free(wt);
	iproc_model_ctx_unref(ctx);
}

double iproc_sloglik_value(iproc_sloglik * sll)
{
	if (!sll)
		return 0.0;

	return (-sll->nsend) * sll->f;
}

/*
 *         g = [ sum{gamma[t,i]} * xbar[0,i]
 *               + ( X[0,i])^T * sum{dP[t,i]}
 *               + sum{dxbar[t,i]} ]
 *             -
 *             [ (X[0,i])^T n[i] + sum{dx[t,i,j]} ]
 */
static void
acc_grad_nocache(struct vector *dst_vector, double scale, iproc_sloglik * sll)
{
	if (!sll)
		return;

	iproc_group_model *group =
	    iproc_model_send_group(sll->model, sll->isend);

	// sum{gamma[t,i]} * xbar[0,i]
	vector_axpy(scale * sll->gamma, group->xbar0, dst_vector);

	// (X[0,i])^T * sum{dP[t,i]}
	design_muls0(scale, TRANS_TRANS,
		     sll->model->design, sll->isend, sll->dp, 1.0, dst_vector);

	// sum{dxbar[t,i]}
	svector_axpy(scale, sll->dxbar, dst_vector);

	// - (X[0,i])^T n[i]
	design_muls0(-scale / sll->nsend, TRANS_TRANS,
		     sll->model->design, sll->isend, sll->nrecv,
		     1.0, dst_vector);

	// -sum{dx[t,i,j]}
	svector_axpy(-scale, sll->dxobs, dst_vector);
}

struct vector *iproc_sloglik_grad(iproc_sloglik * sll)
{
	if (!sll)
		return NULL;

	if (!sll->grad_cached) {
		vector_fill(sll->grad, 0.0);
		acc_grad_nocache(sll->grad, (-sll->nsend) * 1.0, sll);
		sll->grad_cached = true;
	}
	return sll->grad;
}
