#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#include "design.h"

static void iproc_design_clear_svectors(struct darray *svectors)
{
	if (svectors) {
		ssize_t i, n = darray_size(svectors);
		for (i = 0; i < n; i++) {
			struct svector *x =
			    *(struct svector **)darray_at(svectors,
							  i);
			svector_free(x);
		}
		darray_resize(svectors, 0);
	}
}

static void iproc_design_svectors_deinit(struct darray *svectors)
{
	if (svectors) {
		iproc_design_clear_svectors(svectors);
		darray_deinit(svectors);
	}
}

static void iproc_design_ctx_free_dealloc(iproc_design_ctx * ctx)
{
	if (ctx) {
		assert(darray_size(&ctx->dxs) == 0);
		darray_deinit(&ctx->dxs);
		free(ctx);
	}
}

static void iproc_design_ctxs_deinit(struct darray *ctxs)
{
	if (ctxs) {
		ssize_t i, n = darray_size(ctxs);
		for (i = 0; i < n; i++) {
			iproc_design_ctx *ctx =
			    *(iproc_design_ctx **) darray_at(ctxs,
							     i);
			iproc_design_ctx_free_dealloc(ctx);
		}
		darray_deinit(ctxs);
	}
}

static void iproc_design_vars_deinit(struct darray *vars)
{
	if (vars) {
		ssize_t i, n = darray_size(vars);
		for (i = 0; i < n; i++) {
			iproc_design_var *var =
			    *(iproc_design_var **) darray_at(vars,
							     i);
			if (var->free)
				var->free(var);
		}
		darray_deinit(vars);
	}
}

void design_deinit(struct design *design)
{
	assert(design);
	iproc_design_svectors_deinit(&design->svectors);
	iproc_design_ctxs_deinit(&design->ctxs);
	iproc_design_vars_deinit(&design->vars);
	actors_free(design->receivers);
	actors_free(design->senders);
}

void
iproc_design_var_init(iproc_design_var * var,
		      ssize_t dim,
		      void (*get_dxs) (iproc_design_var *,
				       iproc_design_ctx * ctx,
				       ssize_t),
		      void (*free) (iproc_design_var *))
{
	assert(var);
	assert(dim >= 0);

	var->dim = dim;
	var->get_dxs = get_dxs;
	var->free = free;
}

bool design_init(struct design *design, struct actors *senders,
		 struct actors *receivers, bool has_reffects, bool has_loops,
		 const struct vector *intervals)
{
	assert(design);
	assert(senders);
	assert(receivers);
	assert(intervals);

#ifndef NDEBUG
	ssize_t i, nintervals = vector_dim(intervals);
	for (i = 1; i < nintervals; i++) {
		assert(vector_get(intervals, i-1) < vector_get(intervals, i));
	}
#endif
	     
	ssize_t nreceivers = actors_size(receivers);
	ssize_t p = actors_dim(senders);
	ssize_t q = actors_dim(receivers);
	
	if (!darray_init(&design->design_dyad_vars, sizeof(struct design_dyad_var)))
		goto fail_design_dyad_vars;
	
	if (!darray_init(&design->vars, sizeof(iproc_design_var *)))
		goto fail_vars;
	
	if (!darray_init(&design->ctxs, sizeof(iproc_design_ctx *)))
		goto fail_ctxs;
	
	if (!darray_init(&design->svectors, sizeof(struct svector *)))
		goto fail_svectors;
	
	if (!refcount_init(&design->refcount))
		goto fail_refcount;
	
	design->senders = actors_ref(senders);
	design->receivers = actors_ref(receivers);
	design->intervals = intervals;
	design->has_reffects = has_reffects;
	design->ireffects = 0;
	design->nreffects = has_reffects ? nreceivers : 0;
	design->istatic = design->ireffects + design->nreffects;
	design->nstatic = p * q;
	design->idynamic = design->istatic + design->nstatic;
	design->ndynamic = 0;
	design->dim = design->idynamic + design->ndynamic;
	design->has_loops = has_loops;
	return true;
	
fail_refcount:
	darray_deinit(&design->svectors);
fail_svectors:
	darray_deinit(&design->ctxs);
fail_ctxs:
	darray_deinit(&design->vars);
fail_vars:
	darray_deinit(&design->design_dyad_vars);
fail_design_dyad_vars:
	return false;
}

struct design *design_alloc(struct actors *senders, struct actors *receivers,
			    bool has_reffects, bool has_loops,
			    const struct vector *intervals)
{
	struct design *design = malloc(sizeof(*design));

	if (design) {
		if (design_init(design, senders, receivers, has_reffects, has_loops,
				intervals))
			return design;

		free(design);
	}

	return NULL;
}

struct design *design_ref(struct design *design)
{
	assert(design);
	refcount_get(&design->refcount);
	return design;
}

void design_free(struct design * design)
{
	if (design && refcount_put(&design->refcount, NULL)) {
		design_deinit(design);
		free(design);
	}
}

ssize_t design_dim(const struct design * design)
{
	assert(design);
	return design->dim;
}

void iproc_design_append(struct design * design, iproc_design_var * var)
{
	assert(design);
	assert(var);
	assert(var->dim >= 0);
	assert(var->get_dxs);

	// the old svectors are invalid since the dimension has changed
	iproc_design_clear_svectors(&design->svectors);
	darray_push_back(&design->vars, &var);
	design->ndynamic += var->dim;
	design->dim += var->dim;
}

static void
design_mul0_reffects(double alpha,
			   enum trans_op trans,
			   const struct design * design,
			   const struct vector *x, struct vector *y)
{
	if (!design->has_reffects)
		return;

	ssize_t off = design->ireffects;
	ssize_t dim = design->nreffects;

	if (trans == TRANS_NOTRANS) {
		struct vector xsub;
		vector_init_slice(&xsub, x, off, dim);
		vector_axpy(alpha, &xsub, y);
	} else {
		struct vector ysub;
		vector_init_slice(&ysub, y, off, dim);
		vector_axpy(alpha, x, &ysub);
	}
}

static void
design_muls0_reffects(double alpha,
			    enum trans_op trans,
			    const struct design * design,
			    const struct svector *x, struct vector *y)
{
	if (!design->has_reffects)
		return;

	ssize_t off = design->ireffects;
	ssize_t dim = design->nreffects;
	ssize_t end = off + dim;

	if (trans == TRANS_NOTRANS) {
		struct svector_iter itx;
		ssize_t i;
		double x_i;

		svector_iter_init(x, &itx);
		while (svector_iter_advance(x, &itx)) {
			i = svector_iter_current_index(x, &itx);
			if (off <= i && i < end) {
				x_i = *svector_iter_current(x, &itx);
				*vector_at(y, i) += alpha * x_i;
			}
		}
		svector_iter_deinit(x, &itx);
	} else {
		struct vector ysub;
		vector_init_slice(&ysub, y, off, dim);
		svector_axpy(alpha, x, &ysub);
	}
}

static void
design_mul0_static(double alpha,
			 enum trans_op trans,
			 const struct design * design,
			 ssize_t isend, const struct vector *x,
			 struct vector *y)
{
	if (design->nstatic == 0)
		return;

	const struct actors *senders = design_senders(design);
	const struct actors *receivers = design_receivers(design);
	ssize_t p = actors_dim(senders);
	ssize_t q = actors_dim(receivers);
	ssize_t ix_begin = design->istatic;
	ssize_t nstatic = design->nstatic;
	const struct vector *s = actors_traits(senders, isend);
	struct vector *z = vector_alloc(q);

	if (trans == TRANS_NOTRANS) {
		struct vector xsub;
		vector_init_slice(&xsub, x, ix_begin, nstatic);

		/* z := alpha t(x) s */
		struct matrix xmat;
		matrix_init_view_vector(&xmat, &xsub, p, q);

		matrix_mul(alpha, TRANS_TRANS, &xmat, s, 0.0, z);

		/* y := y + R z */
		actors_mul(1.0, TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_mul(alpha, TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		struct vector ysub;
		vector_init_slice(&ysub, y, ix_begin, nstatic);

		struct matrix ymat;
		matrix_init_view_vector(&ymat, &ysub, p, q);

		struct matrix smat;
		matrix_init_view_vector(&smat, s, p, 1);

		struct matrix zmat;
		matrix_init_view_vector(&zmat, z, 1, q);

		matrix_matmul(1.0, TRANS_NOTRANS, &smat, &zmat, 1.0, &ymat);
	}

	vector_free(z);
}

static void
design_muls0_static(double alpha,
			  enum trans_op trans,
			  const struct design * design,
			  ssize_t isend, const struct svector *x,
			  struct vector *y)
{
	if (design->nstatic == 0)
		return;

	const struct actors *senders = design_senders(design);
	const struct actors *receivers = design_receivers(design);
	ssize_t p = actors_dim(senders);
	ssize_t q = actors_dim(receivers);
	ssize_t ix_begin = design->istatic;
	ssize_t nstatic = design->nstatic;
	ssize_t ix_end = ix_begin + nstatic;
	const struct vector *s = actors_traits(senders, isend);
	struct vector *z = vector_alloc(q);

	if (trans == TRANS_NOTRANS) {
		/* z := alpha t(x) s 
		 *
		 * z[j] = alpha * { \sum_i (x[i,j] * s[i]) }
		 */

		struct svector_iter itx;
		double x_ij, s_i;
		ssize_t ix, i, j;
		imaxdiv_t ij;

		vector_fill(z, 0.0);

		svector_iter_init(x, &itx);
		while (svector_iter_advance(x, &itx)) {
			ix = svector_iter_current_index(x, &itx);
			if (ix < ix_begin || ix >= ix_end)
				continue;

			ij = imaxdiv(ix - ix_begin, p);
			i = ij.rem;	/* ix % p */
			j = ij.quot;	/* ix / p */
			x_ij = *svector_iter_current(x, &itx);
			s_i = *vector_at(s, i);
			*vector_at(z, j) += x_ij * s_i;
		}
		svector_iter_deinit(x, &itx);

		vector_scale(z, alpha);

		/* y := y + R z */
		actors_mul(1.0, TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_muls(alpha, TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		struct vector ysub;
		vector_init_slice(&ysub, y, ix_begin, nstatic);

		struct matrix ymat;
		matrix_init_view_vector(&ymat, &ysub, p, q);

		struct matrix smat;
		matrix_init_view_vector(&smat, s, p, 1);

		struct matrix zmat;
		matrix_init_view_vector(&zmat, z, 1, q);

		matrix_matmul(1.0, TRANS_NOTRANS, &smat, &zmat, 1.0, &ymat);
	}

	vector_free(z);
}

void
design_mul0(double alpha,
		  enum trans_op trans,
		  const struct design * design,
		  ssize_t isend,
		  const struct vector *x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < design_nsender(design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_dim(design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_mul0_reffects(alpha, trans, design, x, y);
	design_mul0_static(alpha, trans, design, isend, x, y);
}

void
design_muls0(double alpha,
		   enum trans_op trans,
		   const struct design * design,
		   ssize_t isend,
		   const struct svector *x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < design_nsender(design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_dim(design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_muls0_reffects(alpha, trans, design, x, y);
	design_muls0_static(alpha, trans, design, isend, x, y);
}

ssize_t design_nsender(const struct design *design)
{
	assert(design);
	const struct actors *senders = design_senders(design);
	return actors_size(senders);
}

ssize_t design_nreceiver(const struct design *design)
{
	assert(design);
	const struct actors *receivers = design_receivers(design);
	return actors_size(receivers);
}

struct actors *design_senders(const struct design *design)
{
	assert(design);
	return design->senders;
}

struct actors *design_receivers(const struct design *design)
{
	assert(design);
	return design->receivers;
}

bool design_has_reffects(const struct design *design)
{
	assert(design);
	return design->has_reffects;
}

bool design_has_loops(const struct design *design)
{
	assert(design);
	return design->has_loops;
}

ssize_t design_add_dyad_var(struct design *design, const struct dyad_var *var)
{
	assert(design);
	assert(var);
	assert(var->dim >= 0);
	assert(var->get_jrecv_dxs);
	       
	struct design_dyad_var *design_var;
	
	if ((design_var = darray_push_back(&design->design_dyad_vars, NULL))) {
		design_var->index = design->idynamic + design->ndynamic;
		design_var->var = var;
		design->ndynamic += design_var->var->dim;
		design->dim += design_var->var->dim;
		return design_var->index;
	}

	return -1;
}
