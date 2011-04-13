#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include "memory.h"
#include "sloglik.h"

static void iproc_sloglik_free(iproc_sloglik * sll)
{
	if (sll) {
		iproc_svector_unref(sll->dxbar);
		iproc_svector_unref(sll->dp);
		iproc_svector_unref(sll->dxobs);
		iproc_svector_unref(sll->nrecv);
		vector_free(sll->grad);
		iproc_model_unref(sll->model);
		iproc_free(sll);
	}
}

iproc_sloglik *iproc_sloglik_new(iproc_model * model, int64_t isend)
{
	iproc_sloglik *sll = iproc_malloc(sizeof(*sll));
	if (!sll)
		return NULL;

	int64_t n = iproc_model_nreceiver(model);
	int64_t p = iproc_model_dim(model);

	sll->model = iproc_model_ref(model);
	sll->isend = isend;

	sll->f = 0.0;
	sll->grad = vector_new(p);
	sll->grad_cached = false;

	sll->nsend = 0;
	sll->nrecv = iproc_svector_new(n);
	sll->dxobs = iproc_svector_new(p);

	sll->gamma = 0.0;
	sll->dp = iproc_svector_new(n);
	sll->dxbar = iproc_svector_new(p);

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

void
iproc_sloglik_insert(iproc_sloglik * sll,
		     iproc_history * history, int64_t jrecv)
{
	iproc_sloglik_insertm(sll, history, &jrecv, 1);
}

void
iproc_sloglik_insertm(iproc_sloglik * sll,
		      iproc_history * history, int64_t * jrecv, int64_t n)
{
	int64_t nreceiver = iproc_model_nreceiver(sll->model);
	iproc_model_ctx *ctx =
	    iproc_model_ctx_new(sll->model, sll->isend, history);
	iproc_svector *wt = iproc_svector_new(nreceiver);
	double ntot = sll->nsend + n;
	double scale1 = n / ntot;
	double scale0 = 1 - scale1;
	double lpbar = 0.0;
	int64_t i;

	sll->grad_cached = false;

	for (i = 0; i < n; i++) {
		assert(jrecv[i] >= 0);
		assert(jrecv[i] < nreceiver);

		double lp = iproc_model_ctx_logprob(ctx, jrecv[i]);
		lpbar += (lp - lpbar) / (i + 1);

		iproc_svector_inc(wt, jrecv[i], 1.0);

		// update number of receives
		iproc_svector_inc(sll->nrecv, jrecv[i], 1.0);
	}

	// update log likelihood
	sll->f += scale1 * ((-lpbar) - sll->f);

	// update observed variable diffs
	iproc_design_ctx_dmuls(scale1 / n, IPROC_TRANS_TRANS, ctx->design_ctx,
			       wt, scale0, sll->dxobs);

	sll->gamma += scale1 * (ctx->gamma - sll->gamma);

	iproc_svector_scale(sll->dp, scale0);
	iproc_svector_sacc(sll->dp, scale1, ctx->dp);

	iproc_svector_scale(sll->dxbar, scale0);
	iproc_svector_sacc(sll->dxbar, scale1, ctx->dxbar);

	// update number of sends
	sll->nsend += n;

	iproc_svector_unref(wt);
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
	iproc_design_muls0(scale, IPROC_TRANS_TRANS,
			   sll->model->design, sll->isend, sll->dp,
			   1.0, dst_vector);

	// sum{dxbar[t,i]}
	iproc_vector_sacc(dst_vector, scale, sll->dxbar);

	// - (X[0,i])^T n[i]
	iproc_design_muls0(-scale / sll->nsend, IPROC_TRANS_TRANS,
			   sll->model->design, sll->isend, sll->nrecv,
			   1.0, dst_vector);

	// -sum{dx[t,i,j]}
	iproc_vector_sacc(dst_vector, -scale, sll->dxobs);
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
