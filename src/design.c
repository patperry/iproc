#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#include "design.h"

static void iproc_design_clear_svectors(struct darray *svectors)
{
	if (svectors) {
		int64_t i, n = darray_size(svectors);
		for (i = 0; i < n; i++) {
			struct svector *x =
			    *(struct svector **) darray_at(svectors,
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
		int64_t i, n = darray_size(ctxs);
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
		int64_t i, n = darray_size(vars);
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

static void iproc_design_free(iproc_design * design)
{

	if (design) {
		iproc_design_svectors_deinit(&design->svectors);
		iproc_design_ctxs_deinit(&design->ctxs);
		iproc_design_vars_deinit(&design->vars);
		actors_free(design->receivers);
		actors_free(design->senders);
		free(design);
	}
}

void
iproc_design_var_init(iproc_design_var * var,
		      int64_t dim,
		      void (*get_dxs) (iproc_design_var *,
				       iproc_design_ctx * ctx,
				       int64_t),
		      void (*free) (iproc_design_var *))
{
	assert(var);
	assert(dim >= 0);

	var->dim = dim;
	var->get_dxs = get_dxs;
	var->free = free;
}

iproc_design *iproc_design_new(struct actors *senders,
			       struct actors *receivers, bool has_reffects)
{
	assert(senders);
	assert(receivers);

	iproc_design *design = calloc(1, sizeof(*design));

	if (!design)
		return NULL;

	int64_t nreceivers = actors_size(receivers);
	int64_t p = actors_dim(senders);
	int64_t q = actors_dim(receivers);

	design->senders = actors_ref(senders);
	design->receivers = actors_ref(receivers);
	design->has_reffects = has_reffects;
	design->ireffects = 0;
	design->nreffects = has_reffects ? nreceivers : 0;
	design->istatic = design->ireffects + design->nreffects;
	design->nstatic = p * q;
	design->idynamic = design->istatic + design->nstatic;
	design->ndynamic = 0;
	design->dim = design->idynamic + design->ndynamic;
	refcount_init(&design->refcount);

	if (!(darray_init(&design->vars, sizeof(iproc_design_var *))
	      && darray_init(&design->ctxs, sizeof(iproc_design_ctx *))
	      && darray_init(&design->svectors, sizeof(struct svector *)))) {
		iproc_design_free(design);
		design = NULL;
	}

	return design;
}

iproc_design *iproc_design_ref(iproc_design * design)
{
	if (design) {
		refcount_get(&design->refcount);
	}
	return design;
}

static void iproc_design_release(struct refcount *refcount)
{
	iproc_design *design = container_of(refcount, iproc_design, refcount);
	iproc_design_free(design);
}

void iproc_design_unref(iproc_design * design)
{
	if (!design)
		return;

	refcount_put(&design->refcount, iproc_design_release);
}

int64_t iproc_design_dim(const iproc_design * design)
{
	assert(design);
	return design->dim;
}

void iproc_design_append(iproc_design * design, iproc_design_var * var)
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
iproc_design_mul0_reffects(double alpha,
			   enum trans_op trans,
			   const iproc_design * design,
			   const struct vector *x, struct vector *y)
{
	if (!design->has_reffects)
		return;

	int64_t off = design->ireffects;
	int64_t dim = design->nreffects;

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
iproc_design_muls0_reffects(double alpha,
			    enum trans_op trans,
			    const iproc_design * design,
			    const struct svector * x, struct vector *y)
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
iproc_design_mul0_static(double alpha,
			 enum trans_op trans,
			 const iproc_design * design,
			 int64_t isend, const struct vector *x,
			 struct vector *y)
{
	if (design->nstatic == 0)
		return;

	const struct actors *senders = iproc_design_senders(design);
	const struct actors *receivers = iproc_design_receivers(design);
	int64_t p = actors_dim(senders);
	int64_t q = actors_dim(receivers);
	int64_t ix_begin = design->istatic;
	int64_t nstatic = design->nstatic;
	const struct vector *s = actors_traits(senders, isend);
	struct vector *z = vector_alloc(q);

	if (trans == TRANS_NOTRANS) {
		struct vector xsub;
		vector_init_slice(&xsub, x, ix_begin, nstatic);

		/* z := alpha t(x) s */
		iproc_matrix_view xmat =
		    iproc_matrix_view_vector(&xsub, p, q);
		matrix_mul(alpha, TRANS_TRANS, &xmat.matrix, s, 0.0,
				 z);

		/* y := y + R z */
		actors_mul(1.0, TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_mul(alpha, TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		struct vector ysub;
		vector_init_slice(&ysub, y, ix_begin, nstatic);
		iproc_matrix_view ymat =
		    iproc_matrix_view_vector(&ysub, p, q);
		iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
		iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
		matrix_matmul(1.0, TRANS_NOTRANS, &smat.matrix,
				    &zmat.matrix, 1.0, &ymat.matrix);
	}

	vector_free(z);
}

static void
iproc_design_muls0_static(double alpha,
			  enum trans_op trans,
			  const iproc_design * design,
			  int64_t isend, const struct svector * x,
			  struct vector *y)
{
	if (design->nstatic == 0)
		return;

	const struct actors *senders = iproc_design_senders(design);
	const struct actors *receivers = iproc_design_receivers(design);
	int64_t p = actors_dim(senders);
	int64_t q = actors_dim(receivers);
	int64_t ix_begin = design->istatic;
	int64_t nstatic = design->nstatic;
	int64_t ix_end = ix_begin + nstatic;
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
		iproc_matrix_view ymat =
		    iproc_matrix_view_vector(&ysub, p, q);
		iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
		iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
		matrix_matmul(1.0, TRANS_NOTRANS, &smat.matrix,
				    &zmat.matrix, 1.0, &ymat.matrix);
	}

	vector_free(z);
}

void
iproc_design_mul0(double alpha,
		  enum trans_op trans,
		  const iproc_design * design,
		  int64_t isend,
		  const struct vector *x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < iproc_design_nsender(design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == iproc_design_dim(design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == iproc_design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	iproc_design_mul0_reffects(alpha, trans, design, x, y);
	iproc_design_mul0_static(alpha, trans, design, isend, x, y);
}

void
iproc_design_muls0(double alpha,
		   enum trans_op trans,
		   const iproc_design * design,
		   int64_t isend,
		   const struct svector * x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < iproc_design_nsender(design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == iproc_design_dim(design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == iproc_design_nreceiver(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == iproc_design_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	iproc_design_muls0_reffects(alpha, trans, design, x, y);
	iproc_design_muls0_static(alpha, trans, design, isend, x, y);
}

/////// make these static ?

int64_t iproc_design_nsender(const iproc_design * design)
{
	assert(design);
	const struct actors *senders = iproc_design_senders(design);
	return actors_size(senders);
}

int64_t iproc_design_nreceiver(const iproc_design * design)
{
	assert(design);
	const struct actors *receivers = iproc_design_receivers(design);
	return actors_size(receivers);
}

struct actors *iproc_design_senders(const iproc_design * design)
{
	assert(design);
	return design->senders;
}

struct actors *iproc_design_receivers(const iproc_design * design)
{
	assert(design);
	return design->receivers;
}
