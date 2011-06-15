#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#include "design.h"

void design_deinit(struct design *design)
{
	assert(design);
	refcount_deinit(&design->refcount);
	vector_deinit(&design->intervals);
	actors_free(design->receivers);
	actors_free(design->senders);

	ssize_t i, n = array_count(&design->vars);
	struct design_var *v;
	for (i = 0; i < n; i++) {
		v = array_item(&design->vars, i);
		if (v->type->deinit) {
			v->type->deinit(v);
		}
	}
	array_deinit(&design->vars);
}

void design_init(struct design *design, struct actors *senders,
		 struct actors *receivers, const struct vector *intervals)
{
	assert(design);
	assert(senders);
	assert(receivers);
	assert(intervals);

#ifndef NDEBUG
	ssize_t i, nintervals = vector_dim(intervals);
	for (i = 1; i < nintervals; i++) {
		assert(vector_item(intervals, i - 1) <
		       vector_item(intervals, i));
	}
#endif

	ssize_t p = actors_dim(senders);
	ssize_t q = actors_dim(receivers);

	array_init(&design->vars, sizeof(struct design_var));

	design->senders = actors_ref(senders);
	design->receivers = actors_ref(receivers);

	vector_init_copy(&design->intervals, intervals);
	refcount_init(&design->refcount);

	design->reffects = false;
	design->ireffects = 0;
	design->istatic = 0;
	design->nstatic = p * q;
	design->idynamic = design->istatic + design->nstatic;
	design->ndynamic = 0;
	design->dim = design->idynamic + design->ndynamic;
	design->loops = false;
}

struct design *design_alloc(struct actors *senders, struct actors *receivers,
			    const struct vector *intervals)
{
	struct design *design = xcalloc(1, sizeof(*design));
	design_init(design, senders, receivers, intervals);
	return design;
}

struct design *design_ref(struct design *design)
{
	assert(design);
	refcount_get(&design->refcount);
	return design;
}

void design_free(struct design *design)
{
	if (design && refcount_put(&design->refcount, NULL)) {
		refcount_get(&design->refcount);
		design_deinit(design);
		xfree(design);
	}
}

static void
design_mul0_reffects(double alpha,
		     enum trans_op trans,
		     const struct design *design,
		     const struct vector *x, struct vector *y)
{
	if (!design->reffects)
		return;

	ssize_t off = design->ireffects;
	ssize_t dim = design_receiver_count(design);

	if (trans == TRANS_NOTRANS) {
		struct vector xsub = vector_slice(x, off, dim);
		vector_axpy(alpha, &xsub, y);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		vector_axpy(alpha, x, &ysub);
	}
}

static void
design_muls0_reffects(double alpha,
		      enum trans_op trans,
		      const struct design *design,
		      const struct svector *x, struct vector *y)
{
	if (!design->reffects)
		return;

	ssize_t off = design->ireffects;
	ssize_t dim = design_receiver_count(design);
	ssize_t end = off + dim;

	if (trans == TRANS_NOTRANS) {
		struct svector_iter itx;
		ssize_t i;
		double x_i;

		SVECTOR_FOREACH(itx, x) {
			i = SVECTOR_IDX(itx);
			if (off <= i && i < end) {
				x_i = SVECTOR_VAL(itx);
				*vector_item_ptr(y, i) += alpha * x_i;
			}
		}
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		svector_axpy(alpha, x, &ysub);
	}
}

static void
design_mul0_static(double alpha,
		   enum trans_op trans,
		   const struct design *design,
		   ssize_t isend, const struct vector *x, struct vector *y)
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
		struct vector xsub = vector_slice(x, ix_begin, nstatic);

		/* z := alpha t(x) s */
		struct matrix xmat = matrix_make(&xsub, p, q);

		matrix_mul(alpha, TRANS_TRANS, &xmat, s, 0.0, z);

		/* y := y + R z */
		actors_mul(1.0, TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_mul(alpha, TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		struct vector ysub = vector_slice(y, ix_begin, nstatic);
		struct matrix ymat = matrix_make(&ysub, p, q);
		struct matrix smat = matrix_make(s, p, 1);
		struct matrix zmat = matrix_make(z, 1, q);

		matrix_matmul(1.0, TRANS_NOTRANS, &smat, &zmat, 1.0, &ymat);
	}

	vector_free(z);
}

static void
design_muls0_static(double alpha,
		    enum trans_op trans,
		    const struct design *design,
		    ssize_t isend, const struct svector *x, struct vector *y)
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

		SVECTOR_FOREACH(itx, x) {
			ix = SVECTOR_IDX(itx);
			if (ix < ix_begin || ix >= ix_end)
				continue;

			ij = imaxdiv(ix - ix_begin, p);
			i = ij.rem;	/* ix % p */
			j = ij.quot;	/* ix / p */
			x_ij = SVECTOR_VAL(itx);
			s_i = *vector_item_ptr(s, i);
			*vector_item_ptr(z, j) += x_ij * s_i;
		}

		vector_scale(z, alpha);

		/* y := y + R z */
		actors_mul(1.0, TRANS_NOTRANS, receivers, z, 1.0, y);
	} else {
		/* z := alpha t(R) x */
		actors_muls(alpha, TRANS_TRANS, receivers, x, 0.0, z);

		/* y := y + s \otimes z */
		struct vector ysub = vector_slice(y, ix_begin, nstatic);
		struct matrix ymat = matrix_make(&ysub, p, q);
		struct matrix smat = matrix_make(s, p, 1);
		struct matrix zmat = matrix_make(z, 1, q);

		matrix_matmul(1.0, TRANS_NOTRANS, &smat, &zmat, 1.0, &ymat);
	}

	vector_free(z);
}

void
design_mul0(double alpha,
	    enum trans_op trans,
	    const struct design *design,
	    ssize_t isend,
	    const struct vector *x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < design_sender_count(design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || vector_dim(x) == design_dim(design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_receiver_count(design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_receiver_count(design));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == design_dim(design));

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
	     const struct design *design,
	     ssize_t isend,
	     const struct svector *x, double beta, struct vector *y)
{
	assert(design);
	assert(isend >= 0);
	assert(isend < design_sender_count(design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || svector_dim(x) == design_dim(design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_receiver_count(design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_receiver_count(design));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == design_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_muls0_reffects(alpha, trans, design, x, y);
	design_muls0_static(alpha, trans, design, isend, x, y);
}

void design_set_loops(struct design *design, bool loops)
{
	assert(design);
	design->loops = loops;
}

bool design_reffects(const struct design *design)
{
	assert(design);
	return design->reffects;
}

void design_set_reffects(struct design *design, bool reffects)
{
	assert(design);

	if (design->reffects == reffects)
		return;

	ssize_t nrecv = design_receiver_count(design);
	ssize_t delta = reffects ? nrecv : -nrecv;

	ssize_t i, n = array_count(&design->vars);
	struct design_var *var;
	for (i = 0; i < n; i++) {
		var = array_item(&design->vars, i);
		var->index += delta;
	}

	design->istatic += delta;
	design->idynamic += delta;
	design->dim += delta;
	design->reffects = reffects;
}

void design_add_var(struct design *design, const struct var_type *type)
{
	assert(design);
	assert(type);
	assert(type->init);

	struct design_var *var = array_add(&design->vars, NULL);

	type->init(var, design);
	assert(var->dim >= 0);
	var->index = design->idynamic + design->ndynamic;
	var->type = type;
	design->ndynamic += var->dim;
	design->dim += var->dim;
}

ssize_t design_traits_index(const struct design *design)
{
	assert(design);
	return design->istatic;
}

ssize_t design_reffects_index(const struct design *design)
{
	assert(design);

	if (design_reffects(design))
		return design->ireffects;
	return -1;
}

static bool design_var_equals(const void *p1, const void *p2)
{
	const struct design_var *v1 = p1, *v2 = p2;
	return v1->type == v2->type;
}

ssize_t design_var_index(const struct design *design,
			 const struct var_type *type)
{
	assert(design);
	assert(type);

	struct design_var *key = container_of(&type, struct design_var, type);
	const struct design_var *var =
	    array_find(&design->vars, (predicate_fn) design_var_equals, key);

	if (var) {
		return var->index;
	}

	return -1;
}
