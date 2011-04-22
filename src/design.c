#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include "memory.h"
#include "design.h"

static void iproc_design_clear_svectors(struct darray *svectors)
{
	if (svectors) {
		int64_t i, n = darray_size(svectors);
		for (i = 0; i < n; i++) {
			iproc_svector *x =
			    *(iproc_svector **) darray_at(svectors,
							  i);
			iproc_svector_unref(x);
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
		iproc_free(ctx);
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
		iproc_free(design);
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

iproc_design *iproc_design_new(iproc_actors * senders,
			       iproc_actors * receivers, bool has_reffects)
{
	assert(senders);
	assert(receivers);

	iproc_design *design = iproc_calloc(1, sizeof(*design));

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
	      && darray_init(&design->svectors, sizeof(iproc_svector *)))) {
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
			   iproc_trans trans,
			   const iproc_design * design,
			   const struct vector *x, struct vector *y)
{
	if (!design->has_reffects)
		return;

	int64_t off = design->ireffects;
	int64_t dim = design->nreffects;

	if (trans == IPROC_TRANS_NOTRANS) {
		iproc_vector_view xsub = vector_slice(x, off, dim);
		vector_axpy(alpha, &xsub.vector, y);
	} else {
		iproc_vector_view ysub = vector_slice(y, off, dim);
		vector_axpy(alpha, x, &ysub.vector);
	}
}

static void
iproc_design_muls0_reffects(double alpha,
			    iproc_trans trans,
			    const iproc_design * design,
			    const iproc_svector * x, struct vector *y)
{
	if (!design->has_reffects)
		return;

	int64_t off = design->ireffects;
	int64_t dim = design->nreffects;
	int64_t end = off + dim;

	if (trans == IPROC_TRANS_NOTRANS) {
		int64_t nnz = iproc_svector_nnz(x);
		int64_t inz = iproc_svector_find_nz(x, off);

		if (inz < 0)
			inz = ~inz;

		while (inz < nnz) {
			int64_t i = iproc_svector_nz(x, inz);
			if (i >= end)
				break;

			double x_i = iproc_svector_nz_get(x, inz);
			*vector_at(y, i) += alpha * x_i;

			inz++;
		}
	} else {
		iproc_vector_view ysub = vector_slice(y, off, dim);
		iproc_vector_sacc(&ysub.vector, alpha, x);
	}
}

static void
iproc_design_mul0_static(double alpha,
			 iproc_trans trans,
			 const iproc_design * design,
			 int64_t isend, const struct vector *x,
			 struct vector *y)
{
	if (design->nstatic == 0)
		return;

	const iproc_actors *senders = iproc_design_senders(design);
	const iproc_actors *receivers = iproc_design_receivers(design);
	int64_t p = actors_dim(senders);
	int64_t q = actors_dim(receivers);
	int64_t ix_begin = design->istatic;
	int64_t nstatic = design->nstatic;
	const struct vector *s = actors_traits(senders, isend);
	struct vector *z = vector_new(q);

	if (trans == IPROC_TRANS_NOTRANS) {
		iproc_vector_view xsub = vector_slice(x, ix_begin, nstatic);

		/* z := alpha t(x) s */
		iproc_matrix_view xmat =
		    iproc_matrix_view_vector(&xsub.vector, p, q);
		iproc_matrix_mul(alpha, IPROC_TRANS_TRANS, &xmat.matrix, s, 0.0,
				 z);

		/* y := y + R z */
		actors_mul(1.0, IPROC_TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_mul(alpha, IPROC_TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		iproc_vector_view ysub = vector_slice(y, ix_begin, nstatic);
		iproc_matrix_view ymat =
		    iproc_matrix_view_vector(&ysub.vector, p, q);
		iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
		iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
		iproc_matrix_matmul(1.0, IPROC_TRANS_NOTRANS, &smat.matrix,
				    &zmat.matrix, 1.0, &ymat.matrix);
	}

	vector_free(z);
}

static void
iproc_design_muls0_static(double alpha,
			  iproc_trans trans,
			  const iproc_design * design,
			  int64_t isend, const iproc_svector * x,
			  struct vector *y)
{
	if (design->nstatic == 0)
		return;

	const iproc_actors *senders = iproc_design_senders(design);
	const iproc_actors *receivers = iproc_design_receivers(design);
	int64_t p = actors_dim(senders);
	int64_t q = actors_dim(receivers);
	int64_t ix_begin = design->istatic;
	int64_t nstatic = design->nstatic;
	int64_t ix_end = ix_begin + nstatic;
	const struct vector *s = actors_traits(senders, isend);
	struct vector *z = vector_new(q);

	if (trans == IPROC_TRANS_NOTRANS) {
		/* z := alpha t(x) s 
		 *
		 * z[j] = alpha * { \sum_i (x[i,j] * s[i]) }
		 */
		vector_fill(z, 0.0);
		int64_t inz, nnz = iproc_svector_nnz(x);
		for (inz = 0; inz < nnz; inz++) {
			int64_t ix = iproc_svector_nz(x, inz);

			if (ix < ix_begin)
				continue;

			if (ix >= ix_end)
				break;

			imaxdiv_t ij = imaxdiv(ix - ix_begin, p);
			int64_t i = ij.rem;	/* ix % p */
			int64_t j = ij.quot;	/* ix / p */
			double x_ij = iproc_svector_nz_get(x, inz);
			double s_i = *vector_at(s, i);
			*vector_at(z, j) += x_ij * s_i;
		}
		vector_scale(z, alpha);

		/* y := y + R z */
		actors_mul(1.0, IPROC_TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_muls(alpha, IPROC_TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		iproc_vector_view ysub = vector_slice(y, ix_begin, nstatic);
		iproc_matrix_view ymat =
		    iproc_matrix_view_vector(&ysub.vector, p, q);
		iproc_matrix_view smat = iproc_matrix_view_vector(s, p, 1);
		iproc_matrix_view zmat = iproc_matrix_view_vector(z, 1, q);
		iproc_matrix_matmul(1.0, IPROC_TRANS_NOTRANS, &smat.matrix,
				    &zmat.matrix, 1.0, &ymat.matrix);
	}

	vector_free(z);
}

void
iproc_design_mul0(double alpha,
		  iproc_trans trans,
		  const iproc_design * design,
		  int64_t isend,
		  const struct vector *x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < iproc_design_nsender(design));
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_size(x) == iproc_design_dim(design));
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_size(y) == iproc_design_nreceiver(design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_size(x) == iproc_design_nreceiver(design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_size(y) == iproc_design_dim(design));

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
		   iproc_trans trans,
		   const iproc_design * design,
		   int64_t isend,
		   const iproc_svector * x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < iproc_design_nsender(design));
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || iproc_svector_dim(x) == iproc_design_dim(design));
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_size(y) == iproc_design_nreceiver(design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || iproc_svector_dim(x) == iproc_design_nreceiver(design));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_size(y) == iproc_design_dim(design));

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
	const iproc_actors *senders = iproc_design_senders(design);
	return actors_size(senders);
}

int64_t iproc_design_nreceiver(const iproc_design * design)
{
	assert(design);
	const iproc_actors *receivers = iproc_design_receivers(design);
	return actors_size(receivers);
}

iproc_actors *iproc_design_senders(const iproc_design * design)
{
	assert(design);
	return design->senders;
}

iproc_actors *iproc_design_receivers(const iproc_design * design)
{
	assert(design);
	return design->receivers;
}
