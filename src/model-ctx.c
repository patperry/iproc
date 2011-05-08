#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "logsumexp.h"
#include "refcount.h"
#include "model.h"

static void
compute_weight_changes(const struct frame *f, ssize_t isend,
		       bool has_loops,
		       struct vector *coefs,
		       struct vector *log_p0,
		       struct svector *deta, double *gamma, double *log_gamma)
{
	/* compute the changes in weights */
	frame_dmul(1.0, TRANS_NOTRANS, f, isend, coefs, 0.0, deta);
	if (!has_loops) {
		svector_set(deta, isend, -INFINITY);
	}

	/* compute the scale for the weight differences */
	double lwmax = svector_max(deta);
	double logscale = MAX(0.0, lwmax);
	double invscale = exp(-logscale);

	struct logsumexp pos, neg;
	logsumexp_init(&pos);
	logsumexp_init(&neg);

	/* treat the initial sum of weights as the first positive difference */
	logsumexp_insert(&pos, -logscale);

	/* compute the log sums of the positive and negative differences in weights */
	struct svector_iter it;
	svector_iter_init(deta, &it);
	while (svector_iter_advance(deta, &it)) {
		//      for (i = 0; i < nnz; i++) {
		ssize_t jrecv = svector_iter_current_index(deta, &it);
		double lp0 = *vector_at(log_p0, jrecv);
		double dlw = *svector_iter_current(deta, &it);
		double log_abs_dw;

		/* When w > w0:
		 *   w - w0      = w0 * [exp(dlw) - 1];
		 *               = w0 * {exp[dlw - log(scale)] - 1/scale} * scale
		 *
		 *   log(w - w0) = log(w0) + log{exp[(dlw) - log(scale)] - 1/scale} + log(scale)
		 */
		if (dlw >= 0) {
			log_abs_dw = lp0 + log(exp(dlw - logscale) - invscale);
			logsumexp_insert(&pos, log_abs_dw);
		} else {
			log_abs_dw = lp0 + log(invscale - exp(dlw - logscale));
			logsumexp_insert(&neg, log_abs_dw);
		}
	}
	svector_iter_deinit(deta, &it);

	double log_sum_abs_dw_p = logsumexp_value(&pos);
	double log_sum_abs_dw_n = logsumexp_value(&neg);

	if (log_sum_abs_dw_n > log_sum_abs_dw_p) {
		/* The sum of the weights is positive, so this only happens as a
		 * result of numerical errors */
		log_sum_abs_dw_n = log_sum_abs_dw_p;
		printf
		    ("\nWARNING: numerical errors in model-ctx.c (compute_new_logprobs)"
		     "\nPlease report a bug to the authors" "\n\n");
	}

	double log_sum_w = (log_sum_abs_dw_p
			    + log1p(-exp(log_sum_abs_dw_n - log_sum_abs_dw_p))
			    + logscale);
	assert(log_sum_abs_dw_p >= log_sum_abs_dw_n);

	*gamma = exp(-log_sum_w);
	*log_gamma = -log_sum_w;
}

/* Given log_p0, log_gamma, and deta, compute p_active.  The probability for
 * receiver j is active if x[t,i,j] is nonzero or j = i and self-loops are
 * dis-allowed.  We rely on the fact that deta and p_active have the same
 * sparsity pattern.
 */
static void
compute_active_probs(struct vector *log_p0,
		     double log_gamma,
		     struct svector *deta, struct svector *p_active)
{
	ssize_t jrecv;
	double lp0_j, lp_j, deta_j, *p_j;
	struct svector_iter it;

	svector_assign_copy(p_active, deta);

	svector_iter_init(p_active, &it);
	while (svector_iter_advance(p_active, &it)) {
		jrecv = svector_iter_current_index(p_active, &it);
		lp0_j = *vector_at(log_p0, jrecv);
		p_j = svector_iter_current(p_active, &it);
		deta_j = *p_j;
		lp_j = MIN(0.0, log_gamma + lp0_j + deta_j);
		*p_j = exp(lp_j);
	}
	svector_iter_deinit(p_active, &it);
}

static void
compute_prob_diffs(struct vector *p0, double gamma, struct svector *p_active)
{
	ssize_t jrecv;
	double p0_j, dp_j, *p_j;
	struct svector_iter it;

	svector_iter_init(p_active, &it);
	while (svector_iter_advance(p_active, &it)) {
		jrecv = svector_iter_current_index(p_active, &it);
		p0_j = *vector_at(p0, jrecv);
		p_j = svector_iter_current(p_active, &it);
		dp_j = *p_j - gamma * p0_j;
		*p_j = dp_j;

	}
	svector_iter_deinit(p_active, &it);
}

static void iproc_model_ctx_free(iproc_model_ctx * ctx)
{
	if (ctx) {
		iproc_model *model = ctx->model;
		darray_push_back(&model->ctxs, &ctx);
		iproc_model_unref(model);
	}
}

static iproc_model_ctx *iproc_model_ctx_new_alloc(iproc_model * model,
						  const struct frame *f,
						  ssize_t isend)
{
	assert(model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(model));

	iproc_model_ctx *ctx = malloc(sizeof(*ctx));
	if (!ctx)
		return NULL;

	ssize_t nreceiver = iproc_model_nreceiver(model);
	ssize_t dim = iproc_model_dim(model);

	ctx->model = iproc_model_ref(model);
	ctx->frame = NULL;
	ctx->group = NULL;

	ctx->deta = svector_alloc(nreceiver);
	ctx->dp = svector_alloc(nreceiver);
	ctx->dxbar = svector_alloc(dim);

	refcount_init(&ctx->refcount);

	if (!(ctx->deta && ctx->dp && ctx->dxbar)) {
		iproc_model_ctx_free(ctx);
		return NULL;
	}

	iproc_model_ctx_set(ctx, f, isend);
	return ctx;
}

iproc_model_ctx *iproc_model_ctx_new(iproc_model * model,const struct frame *f,
				     ssize_t isend)
{
	assert(model);
	assert(0 <= isend);
	assert(isend < iproc_model_nsender(model));

	iproc_model_ctx *ctx;
	struct darray *ctxs = &model->ctxs;

	if (!darray_empty(ctxs)) {
		ctx = *(iproc_model_ctx **) darray_back(ctxs);
		darray_pop_back(ctxs);
		iproc_model_ref(model);
		refcount_init(&ctx->refcount);
		ctx->frame = NULL;
		iproc_model_ctx_set(ctx, f, isend);
		assert(ctx->model == model);
	} else {
		ctx = iproc_model_ctx_new_alloc(model, f, isend);
	}

	return ctx;

}

void
iproc_model_ctx_set(iproc_model_ctx * ctx, const struct frame *f, ssize_t isend)
{
	assert(ctx);
	assert(ctx->model);
	assert(isend >= 0);
	assert(isend < iproc_model_nsender(ctx->model));

	svector_clear(ctx->deta);
	svector_clear(ctx->dp);
	svector_clear(ctx->dxbar);

	iproc_model *model = ctx->model;
	struct vector *coefs = iproc_model_coefs(model);
	bool has_loops = iproc_model_has_loops(model);
	iproc_group_model *group = iproc_model_send_group(model, isend);

	ctx->frame = f;
	ctx->isend = isend;
	ctx->group = group;

	compute_weight_changes(f, isend, has_loops, coefs, group->log_p0,
			       ctx->deta, &ctx->gamma, &ctx->log_gamma);
	compute_active_probs(group->log_p0, ctx->log_gamma, ctx->deta, ctx->dp);
	frame_dmuls(1.0, TRANS_TRANS, f, isend, ctx->dp, 0.0, ctx->dxbar);
	compute_prob_diffs(group->p0, ctx->gamma, ctx->dp);
}

iproc_model_ctx *iproc_model_ctx_ref(iproc_model_ctx * ctx)
{
	if (ctx) {
		refcount_get(&ctx->refcount);
	}
	return ctx;
}

static void iproc_model_ctx_release(struct refcount *refcount)
{
	iproc_model_ctx *ctx =
	    container_of(refcount, iproc_model_ctx, refcount);
	iproc_model_ctx_free(ctx);
}

void iproc_model_ctx_unref(iproc_model_ctx * ctx)
{
	if (ctx) {
		refcount_put(&ctx->refcount, iproc_model_ctx_release);
	}
}

ssize_t iproc_model_ctx_nreceiver(iproc_model_ctx * ctx)
{
	assert(ctx);
	iproc_model *model = ctx->model;
	ssize_t nreceiver = iproc_model_nreceiver(model);
	return nreceiver;
}

double iproc_model_ctx_prob(iproc_model_ctx * ctx, ssize_t jrecv)
{
	assert(ctx);
	assert(0 <= jrecv);
	assert(jrecv < iproc_model_ctx_nreceiver(ctx));

	double lp = iproc_model_ctx_logprob(ctx, jrecv);
	double p = exp(lp);
	return p;
}

/*
 *         log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 */
double iproc_model_ctx_logprob(iproc_model_ctx * ctx, ssize_t jrecv)
{
	assert(ctx);
	assert(0 <= jrecv);
	assert(jrecv < iproc_model_ctx_nreceiver(ctx));

	/*
	   double gamma = ctx->gamma;
	   double p0 = *vector_at(ctx->group->p0, jrecv);
	   double dp = svector_get(ctx->dp, jrecv);
	   double p = gamma * p0 + dp;
	   p = MAX(0.0, MIN(1.0, p));
	   return log(p);
	 */

	double log_gamma = ctx->log_gamma;
	double log_p0 = *vector_at(ctx->group->log_p0, jrecv);
	double deta = svector_get(ctx->deta, jrecv);
	double log_p = log_gamma + log_p0 + deta;
	return MIN(log_p, 0.0);
}

void iproc_model_ctx_get_probs(iproc_model_ctx * ctx, struct vector *probs)
{
	assert(ctx);
	assert(probs);
	assert(iproc_model_ctx_nreceiver(ctx) == vector_dim(probs));

	iproc_model_ctx_get_logprobs(ctx, probs);
	vector_exp(probs);
}

void
iproc_model_ctx_get_logprobs(iproc_model_ctx * ctx, struct vector *logprobs)
{
	assert(ctx);
	assert(logprobs);
	assert(iproc_model_ctx_nreceiver(ctx) == vector_dim(logprobs));

	ssize_t j, n = iproc_model_ctx_nreceiver(ctx);

	for (j = 0; j < n; j++) {
		double lp = iproc_model_ctx_logprob(ctx, j);
		*vector_at(logprobs, j) = lp;
	}
}
