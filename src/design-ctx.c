#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#include "compare.h"
#include "design.h"

DEFINE_COMPARE_FN(int64_compare, int64_t)

static iproc_svector *iproc_design_var_new_alloc(iproc_design * design)
{
	return svector_new(design->dim);
}

// TODO : make static
static iproc_svector *iproc_design_var_new(iproc_design * design)
{
	assert(design);
	struct darray *svectors = &design->svectors;
	int64_t n = darray_size(svectors);
	iproc_svector *svector;

	if (n > 0) {
		svector = *(iproc_svector **) darray_at(svectors, n - 1);
		darray_resize(svectors, n - 1);
	} else {
		svector = iproc_design_var_new_alloc(design);
	}

	return svector;
}

// TODO : make static
static void
iproc_sdesign_var_free(iproc_design * design, iproc_svector * svector)
{
	assert(design);
	struct darray *svectors = &design->svectors;

	if (!svector)
		return;

	svector_clear(svector);
	darray_push_back(svectors, &svector);
}

static void clear_dxs(iproc_design * design, struct darray *dxs)
{
	if (!dxs)
		return;

	int64_t i, n = darray_size(dxs);
	for (i = 0; i < n; i++) {
		iproc_design_dx *sv = darray_at(dxs, i);
		iproc_sdesign_var_free(design, sv->dx);
	}

	darray_resize(dxs, 0);
}

static void
iproc_design_ctx_set(iproc_design_ctx * ctx, int64_t isend, iproc_history * h)
{
	assert(ctx);
	assert(ctx->design);
	assert(0 <= isend);
	assert(isend < iproc_design_nsender(ctx->design));

	iproc_design *design = ctx->design;

	if (h != ctx->history) {
		iproc_history_ref(h);
		iproc_history_unref(ctx->history);
		ctx->history = h;
	}
	ctx->isend = isend;
	clear_dxs(ctx->design, &ctx->dxs);

	if (design->ndynamic == 0)
		return;

	struct darray *vars = &design->vars;
	int64_t i, n = darray_size(vars);
	int64_t offset = design->idynamic;

	for (i = 0; i < n; i++) {
		iproc_design_var *var = *(iproc_design_var **) darray_at(vars,
									 i);
		assert(var->get_dxs);
		var->get_dxs(var, ctx, offset);
		offset += var->dim;
	}

	assert(offset == design->idynamic + design->ndynamic);
}

static void iproc_design_ctx_free(iproc_design_ctx * ctx)
{
	if (ctx) {
		iproc_history_unref(ctx->history);
		ctx->history = NULL;
		clear_dxs(ctx->design, &ctx->dxs);

		iproc_design *design = ctx->design;
		darray_push_back(&design->ctxs, &ctx);
		iproc_design_unref(design);
	}

}

static iproc_design_ctx *iproc_design_ctx_new_alloc(iproc_design * design,
						    int64_t isend,
						    iproc_history * h)
{
	assert(design);
	assert(0 <= isend);
	assert(isend < iproc_design_nsender(design));

	iproc_design_ctx *ctx = calloc(1, sizeof(*ctx));

	if (!ctx)
		return NULL;

	ctx->design = iproc_design_ref(design);
	ctx->history = NULL;
	ctx->isend = -1;
	refcount_init(&ctx->refcount);

	if (!darray_init(&ctx->dxs, sizeof(iproc_design_dx))) {
		iproc_design_ctx_free(ctx);
		ctx = NULL;
	} else {
		iproc_design_ctx_set(ctx, isend, h);
	}

	return ctx;
}

iproc_design_ctx *iproc_design_ctx_new(iproc_design * design,
				       int64_t isend, iproc_history * h)
{
	assert(design);
	assert(0 <= isend);
	assert(isend < iproc_design_nsender(design));

	iproc_design_ctx *ctx;
	struct darray *ctxs = &design->ctxs;
	int64_t n = darray_size(ctxs);

	if (n > 0) {
		ctx = *(iproc_design_ctx **) darray_at(ctxs, n - 1);
		darray_resize(ctxs, n - 1);
		iproc_design_ref(design);
		refcount_init(&ctx->refcount);
		iproc_design_ctx_set(ctx, isend, h);
		assert(ctx->design == design);
	} else {
		ctx = iproc_design_ctx_new_alloc(design, isend, h);
	}

	return ctx;
}

iproc_design_ctx *iproc_design_ctx_ref(iproc_design_ctx * ctx)
{
	if (ctx) {
		refcount_get(&ctx->refcount);
	}
	return ctx;
}

static void iproc_design_ctx_release(struct refcount *refcount)
{
	iproc_design_ctx *ctx =
	    //(void *)((char *)refcount - offsetof(iproc_design_ctx, refcount));
	    container_of(refcount, iproc_design_ctx, refcount);
	iproc_design_ctx_free(ctx);
}

void iproc_design_ctx_unref(iproc_design_ctx * ctx)
{
	if (!ctx)
		return;

	refcount_put(&ctx->refcount, iproc_design_ctx_release);
}

iproc_svector *iproc_design_ctx_dx(iproc_design_ctx * ctx,
				   int64_t jrecv, bool null_ok)
{
	assert(ctx);
	assert(jrecv >= 0);
	assert(jrecv < iproc_design_nreceiver(ctx->design));

	struct darray *dxs = &ctx->dxs;
	iproc_svector *dx = NULL;

	int64_t i = darray_binary_search(dxs, &jrecv, int64_compare);

	if (i < 0) {
		if (!null_ok) {
			i = ~i;
			dx = iproc_design_var_new(ctx->design);
			iproc_design_dx new_dx = { jrecv, dx };
			darray_insert(dxs, i, &new_dx);
		}
	} else {
		dx = ((iproc_design_dx *) darray_at(dxs, i))->dx;
	}

	return dx;
}

int64_t iproc_design_ctx_nnz(iproc_design_ctx * ctx)
{
	assert(ctx);
	return darray_size(&ctx->dxs);
}

iproc_svector *iproc_design_ctx_nz(iproc_design_ctx * ctx,
				   int64_t inz, int64_t *jrecv)
{
	assert(inz >= 0);
	assert(inz < iproc_design_ctx_nnz(ctx));
	assert(jrecv);

	iproc_design_dx dx = *(iproc_design_dx *) darray_at(&ctx->dxs,
							    inz);
	*jrecv = dx.jrecv;
	return dx.dx;
}

void
iproc_design_ctx_mul(double alpha,
		     iproc_trans trans,
		     iproc_design_ctx * ctx,
		     struct vector *x, double beta, struct vector *y)
{
	assert(ctx);
	assert(ctx->design);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_dim(x) == iproc_design_dim(ctx->design));
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_dim(x) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_dim(ctx->design));

	iproc_design_mul0(alpha, trans, ctx->design, ctx->isend, x, beta, y);

	iproc_svector *diffprod = svector_new(vector_dim(y));
	iproc_design_ctx_dmul(alpha, trans, ctx, x, 0.0, diffprod);
	svector_axpy(1.0, diffprod, y);
	svector_free(diffprod);
}

void
iproc_design_ctx_muls(double alpha,
		      iproc_trans trans,
		      iproc_design_ctx * ctx,
		      iproc_svector * x, double beta, struct vector *y)
{
	assert(ctx);
	assert(ctx->design);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || svector_dim(x) == iproc_design_dim(ctx->design));
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || svector_dim(x) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_dim(ctx->design));

	iproc_design_muls0(alpha, trans, ctx->design, ctx->isend, x, beta, y);

	iproc_svector *diffprod = svector_new(vector_dim(y));
	iproc_design_ctx_dmuls(alpha, trans, ctx, x, 0.0, diffprod);
	svector_axpy(1.0, diffprod, y);
	svector_free(diffprod);
}

void
iproc_design_ctx_dmul(double alpha,
		      iproc_trans trans,
		      const iproc_design_ctx * ctx,
		      const struct vector *x, double beta, iproc_svector * y)
{
	assert(ctx);
	assert(ctx->design);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_dim(x) == iproc_design_dim(ctx->design));
	assert(trans != IPROC_TRANS_NOTRANS
	       || svector_dim(y) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_dim(x) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || svector_dim(y) == iproc_design_dim(ctx->design));

	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}

	if (ctx->history == NULL)
		return;

	iproc_design *design = ctx->design;
	int64_t ndynamic = design->ndynamic;

	if (ndynamic == 0)
		return;

	const struct darray *dxs = &ctx->dxs;

	if (trans == IPROC_TRANS_NOTRANS) {
		int64_t i, n = darray_size(dxs);

		for (i = 0; i < n; i++) {
			iproc_design_dx *sv = darray_at(dxs, i);
			int64_t jrecv = sv->jrecv;
			iproc_svector *dx = sv->dx;
			double dot = svector_dot(dx, x);
			iproc_svector_inc(y, jrecv, alpha * dot);
		}
	} else {
		int64_t i, n = darray_size(dxs);

		for (i = 0; i < n; i++) {
			iproc_design_dx *sv = darray_at(dxs, i);
			int64_t jrecv = sv->jrecv;
			iproc_svector *dx = sv->dx;
			double xjrecv = *vector_at(x, jrecv);

			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			svector_axpys(jscale, dx, y);
		}
	}
}

void
iproc_design_ctx_dmuls(double alpha,
		       iproc_trans trans,
		       const iproc_design_ctx * ctx,
		       const iproc_svector * x, double beta, iproc_svector * y)
{
	assert(ctx);
	assert(ctx->design);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || svector_dim(x) == iproc_design_dim(ctx->design));
	assert(trans != IPROC_TRANS_NOTRANS
	       || svector_dim(y) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || svector_dim(x) == iproc_design_nreceiver(ctx->design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || svector_dim(y) == iproc_design_dim(ctx->design));

	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}

	if (ctx->history == NULL)
		return;

	iproc_design *design = ctx->design;
	int64_t ndynamic = design->ndynamic;

	if (ndynamic == 0)
		return;

	int64_t ix_begin = design->idynamic;
	int64_t ix_end = ix_begin + ndynamic;
	const struct darray *dxs = &ctx->dxs;

	if (trans == IPROC_TRANS_NOTRANS) {
		int64_t i, n = darray_size(dxs);

		for (i = 0; i < n; i++) {
			iproc_design_dx *sv = darray_at(dxs,
							i);
			int64_t jrecv = sv->jrecv;
			iproc_svector *dx = sv->dx;
			int64_t inz, nnz = svector_size(x);
			double dot = 0.0;
			for (inz = 0; inz < nnz; inz++) {
				int64_t ix = iproc_svector_nz(x, inz);
				if (ix < ix_begin || ix >= ix_end)
					continue;

				double xval = iproc_svector_nz_get(x, inz);
				double diffval = svector_get(dx, ix);
				dot += xval * diffval;
			}

			iproc_svector_inc(y, jrecv, alpha * dot);
		}
	} else {
		int64_t i, n = darray_size(dxs);

		for (i = 0; i < n; i++) {
			iproc_design_dx *sv = darray_at(dxs,
							i);
			int64_t jrecv = sv->jrecv;
			iproc_svector *dx = sv->dx;
			double xjrecv = svector_get(x, jrecv);
			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * diff[j] */
			double jscale = alpha * xjrecv;
			svector_axpys(jscale, dx, y);
		}
	}

}
