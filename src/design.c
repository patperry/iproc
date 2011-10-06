#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdlib.h>
#include "vars.h"
#include "design.h"

void design_deinit(struct design *design)
{
	assert(design);
	struct design_var *v;

	ARRAY_FOREACH(v, &design->recv_vars) {
		if (v->type->deinit) {
			v->type->deinit(v);
		}
	}
	array_deinit(&design->recv_vars);

	ARRAY_FOREACH(v, &design->send_vars) {
		if (v->type->deinit) {
			v->type->deinit(v);
		}
	}
	array_deinit(&design->send_vars);

	vector_deinit(&design->intervals);
}

void design_init(struct design *design, struct actors *senders,
		 struct actors *receivers, const struct matrix *traits,
		 const char * const *trait_names,
		 const struct vector *intervals)
{
	assert(design);
	assert(senders);
	assert(receivers);
	assert(traits);
	assert(actors_cohort_count(receivers) == matrix_nrow(traits));
	assert(intervals);

#ifndef NDEBUG
	ssize_t i, nintervals = vector_dim(intervals);
	for (i = 1; i < nintervals; i++) {
		assert(vector_item(intervals, i - 1) <
		       vector_item(intervals, i));
	}
#endif

	ssize_t ps = 0;		// actors_dim(senders);
	ssize_t pr = matrix_ncol(traits);

	design->senders = senders;
	design->receivers = receivers;
	design->traits = traits;
	design->trait_names = trait_names;
	design->loops = false;

	vector_init_copy(&design->intervals, intervals);

	array_init(&design->send_vars, sizeof(struct design_var));
	design->seffects = false;
	design->iseffects = 0;
	design->isstatic = 0;
	design->nsstatic = ps;
	design->isdynamic = design->isstatic + design->nsstatic;
	design->nsdynamic = 0;
	design->sdim = design->isdynamic + design->nsdynamic;

	array_init(&design->recv_vars, sizeof(struct design_var));
	design->reffects = false;
	design->ireffects = 0;
	design->irstatic = 0;
	design->nrstatic = pr;
	design->irdynamic = design->irstatic + design->nrstatic;
	design->nrdynamic = 0;
	design->rdim = design->irdynamic + design->nrdynamic;
}

static void
design_mul0_effects(double alpha,
		    enum blas_trans trans,
		    ssize_t ieffects, ssize_t neffects,
		    const struct vector *x, struct vector *y)
{
	ssize_t off = ieffects;
	ssize_t dim = neffects;

	if (trans == BLAS_NOTRANS) {
		struct vector xsub = vector_slice(x, off, dim);
		vector_axpy(alpha, &xsub, y);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		vector_axpy(alpha, x, &ysub);
	}
}

/*
static void
design_send_mul0_effects(double alpha,
			 enum blas_trans trans,
			 const struct design *design,
			 const struct vector *x, struct vector *y)
{
	if (!design->seffects)
		return;

	ssize_t off = design->iseffects;
	ssize_t dim = design_send_count(design);
	design_mul0_effects(alpha, trans, off, dim, x, y);
}
 */

static void
design_recv_mul0_effects(double alpha,
			 enum blas_trans trans,
			 const struct design *design,
			 const struct vector *x, struct vector *y)
{
	if (!design->reffects)
		return;

	ssize_t off = design->ireffects;
	ssize_t dim = design_recv_count(design);
	design_mul0_effects(alpha, trans, off, dim, x, y);
}

static void
design_muls0_effects(double alpha,
		     enum blas_trans trans,
		     ssize_t ieffects, ssize_t neffects,
		     const struct svector *x, struct vector *y)
{
	ssize_t off = ieffects;
	ssize_t dim = neffects;
	ssize_t end = off + dim;

	if (trans == BLAS_NOTRANS) {
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

/*
static void
design_send_muls0_effects(double alpha,
			  enum blas_trans trans,
			  const struct design *design,
			  const struct svector *x, struct vector *y)
{
	if (!design->seffects)
		return;

	ssize_t off = design->iseffects;
	ssize_t dim = design_send_count(design);
	design_muls0_effects(alpha, trans, off, dim, x, y);
}
*/

static void
design_recv_muls0_effects(double alpha,
			  enum blas_trans trans,
			  const struct design *design,
			  const struct svector *x, struct vector *y)
{
	if (!design->reffects)
		return;

	ssize_t off = design->ireffects;
	ssize_t dim = design_recv_count(design);
	design_muls0_effects(alpha, trans, off, dim, x, y);
}

/*
static void
design_send_mul0_static(double alpha,
			enum blas_trans trans,
			const struct design *design,
			const struct vector *x, struct vector *y)
{
	if (design->nsstatic == 0)
		return;

	const struct actors *senders = design_senders(design);
	ssize_t off = design->isstatic;
	ssize_t dim = design->nsstatic;

	if (trans == BLAS_NOTRANS) {
		struct vector xsub = vector_slice(x, off, dim);
		actors_mul(alpha, BLAS_NOTRANS, senders, &xsub, 1.0, y);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		actors_mul(alpha, TRANS_TRANS, senders, x, 1.0, &ysub);
	}
}
 */

static void
design_recv_mul0_static(double alpha,
			enum blas_trans trans,
			const struct design *design,
			const struct vector *x, struct vector *y)
{
	if (design->nrstatic == 0)
		return;

	const struct actors *receivers = design_receivers(design);
	const struct matrix *traits = design_traits(design);
	ssize_t off = design->irstatic;
	ssize_t dim = design->nrstatic;
	assert(dim == matrix_ncol(traits));
	assert(actors_cohort_count(receivers) == matrix_nrow(traits));

	struct vector z;
	vector_init(&z, actors_cohort_count(receivers));

	if (trans == BLAS_NOTRANS) {
		struct vector xsub = vector_slice(x, off, dim);

		matrix_mul(1.0, trans, traits, &xsub, 0.0, &z);
		actors_mul(alpha, trans, receivers, &z, 1.0, y);

	} else {
		struct vector ysub = vector_slice(y, off, dim);
		actors_mul(1.0, trans, receivers, x, 0.0, &z);
		matrix_mul(alpha, trans, traits, &z, 1.0, &ysub);
	}

	vector_deinit(&z);

	/*
	   const struct actors *senders = design_senders(design);        
	   ssize_t p = actors_dim(senders);
	   ssize_t q = actors_dim(receivers);
	   ssize_t ix_begin = design->irstatic;
	   ssize_t nstatic = design->nrstatic;
	   const struct vector *s = actors_traits(senders, isend);
	   struct vector *z = vector_alloc(q);

	   if (trans == BLAS_NOTRANS) {
	   struct vector xsub = vector_slice(x, ix_begin, nstatic);

	   // z := alpha t(x) s
	   struct matrix xmat = matrix_make(&xsub, p, q);

	   matrix_mul(alpha, TRANS_TRANS, &xmat, s, 0.0, z);

	   // y := y + R z
	   actors_mul(1.0, BLAS_NOTRANS, receivers, z, 1.0, y);
	   } else {
	   // z := alpha t(R) x
	   actors_mul(alpha, TRANS_TRANS, receivers, x, 0.0, z);

	   // y := y + s \otimes z
	   struct vector ysub = vector_slice(y, ix_begin, nstatic);
	   struct matrix ymat = matrix_make(&ysub, p, q);
	   struct matrix smat = matrix_make(s, p, 1);
	   struct matrix zmat = matrix_make(z, 1, q);

	   matrix_matmul(1.0, BLAS_NOTRANS, &smat, &zmat, 1.0, &ymat);
	   }

	   vector_free(z);
	 */
}

/*
static void
design_send_muls0_static(double alpha,
			 enum blas_trans trans,
			 const struct design *design,
			 const struct svector *x, struct vector *y)
{
	if (design->nsstatic == 0)
		return;

	const struct actors *senders = design_senders(design);
	ssize_t off = design->isstatic;
	ssize_t dim = design->nsstatic;

	if (trans == BLAS_NOTRANS) {
		struct svector xsub;
		struct svector_iter it;
		ssize_t i;

		svector_init(&xsub, dim);
		SVECTOR_FOREACH(it, x) {
			i = SVECTOR_IDX(it) - off;
			if (i < 0 || i >= dim)
				continue;
			svector_set_item(&xsub, i, SVECTOR_VAL(it));
		}
		actors_muls(alpha, BLAS_NOTRANS, senders, &xsub, 1.0, y);
		svector_deinit(&xsub);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		actors_muls(alpha, TRANS_TRANS, senders, x, 1.0, &ysub);
	}
}
 */

static void
design_recv_muls0_static(double alpha,
			 enum blas_trans trans,
			 const struct design *design,
			 const struct svector *x, struct vector *y)
{
	if (design->nrstatic == 0)
		return;

	const struct actors *receivers = design_receivers(design);
	const struct matrix *traits = design_traits(design);
	ssize_t off = design->irstatic;
	ssize_t dim = design->nrstatic;
	assert(actors_cohort_count(receivers) == matrix_nrow(traits));
	assert(dim == matrix_ncol(traits));

	if (trans == BLAS_NOTRANS) {
		struct vector z;
		vector_init(&z, actors_cohort_count(receivers));

		struct svector xsub;
		struct svector_iter it;
		ssize_t i;

		svector_init(&xsub, dim);
		SVECTOR_FOREACH(it, x) {
			i = SVECTOR_IDX(it) - off;
			if (i < 0 || i >= dim)
				continue;
			svector_set_item(&xsub, i, SVECTOR_VAL(it));
		}

		matrix_muls(1.0, trans, traits, &xsub, 0.0, &z);
		actors_mul(alpha, trans, receivers, &z, 1.0, y);
		svector_deinit(&xsub);
		vector_deinit(&z);
	} else {
		struct svector z;
		svector_init(&z, actors_cohort_count(receivers));

		struct vector ysub = vector_slice(y, off, dim);
		actors_muls(1.0, trans, receivers, x, 0.0, &z);
		matrix_muls(alpha, trans, traits, &z, 1.0, &ysub);

		svector_deinit(&z);
	}

	/*
	   const struct actors *receivers = design_receivers(design);
	   const struct actors *senders = design_senders(design);       
	   ssize_t p = actors_dim(senders);
	   ssize_t q = actors_dim(receivers);
	   ssize_t off = design->irstatic;
	   ssize_t nstatic = design->nrstatic;
	   ssize_t ix_end = off + nstatic;
	   const struct vector *s = actors_traits(senders, isend);
	   struct vector *z = vector_alloc(q);

	   if (trans == BLAS_NOTRANS) {
	   // z := alpha t(x) s 
	   // z[j] = alpha * { \sum_i (x[i,j] * s[i]) }

	   struct svector_iter itx;
	   double x_ij, s_i;
	   ssize_t ix, i, j;
	   imaxdiv_t ij;

	   vector_fill(z, 0.0);

	   SVECTOR_FOREACH(itx, x) {
	   ix = SVECTOR_IDX(itx);
	   if (ix < off || ix >= ix_end)
	   continue;

	   ij = imaxdiv(ix - off, p);
	   i = (ssize_t)ij.rem; // ix % p
	   j = (ssize_t)ij.quot;        // ix / p
	   x_ij = SVECTOR_VAL(itx);
	   s_i = *vector_item_ptr(s, i);
	   *vector_item_ptr(z, j) += x_ij * s_i;
	   }

	   vector_scale(z, alpha);

	   // y := y + R z
	   actors_mul(1.0, BLAS_NOTRANS, receivers, z, 1.0, y);
	   } else {
	   // z := alpha t(R) x
	   actors_muls(alpha, TRANS_TRANS, receivers, x, 0.0, z);

	   // y := y + s \otimes z
	   struct vector ysub = vector_slice(y, off, nstatic);
	   struct matrix ymat = matrix_make(&ysub, p, q);
	   struct matrix smat = matrix_make(s, p, 1);
	   struct matrix zmat = matrix_make(z, 1, q);

	   matrix_matmul(1.0, BLAS_NOTRANS, &smat, &zmat, 1.0, &ymat);
	   }

	   vector_free(z);
	 */
}

/*
void
design_send_mul0(double alpha,
		 enum blas_trans trans,
		 const struct design *design,
		 const struct vector *x, double beta, struct vector *y)
{
	assert(design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || vector_dim(x) == design_send_dim(design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_send_count(design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(x) == design_send_count(design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_send_dim(design));

	// y := beta y
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_send_mul0_effects(alpha, trans, design, x, y);
	design_send_mul0_static(alpha, trans, design, x, y);
}
*/

void
design_recv_mul0(double alpha,
		 enum blas_trans trans,
		 const struct design *design,
		 const struct vector *x, double beta, struct vector *y)
{
	assert(design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || vector_dim(x) == design_recv_dim(design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_recv_count(design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(x) == design_recv_count(design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_recv_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_recv_mul0_effects(alpha, trans, design, x, y);
	design_recv_mul0_static(alpha, trans, design, x, y);
}

/*
void
design_send_muls0(double alpha,
		  enum blas_trans trans,
		  const struct design *design,
		  const struct svector *x, double beta, struct vector *y)
{
	assert(design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || svector_dim(x) == design_send_dim(design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_send_count(design));
	assert(trans == BLAS_NOTRANS
	       || svector_dim(x) == design_send_count(design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_send_dim(design));

	// y := beta y
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_send_muls0_effects(alpha, trans, design, x, y);
	design_send_muls0_static(alpha, trans, design, x, y);
}
*/

void
design_recv_muls0(double alpha,
		  enum blas_trans trans,
		  const struct design *design,
		  const struct svector *x, double beta, struct vector *y)
{
	assert(design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || svector_dim(x) == design_recv_dim(design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_recv_count(design));
	assert(trans == BLAS_NOTRANS
	       || svector_dim(x) == design_recv_count(design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_recv_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_recv_muls0_effects(alpha, trans, design, x, y);
	design_recv_muls0_static(alpha, trans, design, x, y);
}

void design_set_loops(struct design *design, bool loops)
{
	assert(design);
	design->loops = loops;
}

bool design_send_effects(const struct design *design)
{
	assert(design);
	return design->seffects;
}

bool design_recv_effects(const struct design *design)
{
	assert(design);
	return design->reffects;
}

static void design_set_effects(bool *peffects, bool effects, ssize_t neffects,
			       struct array *vars, ssize_t *istatic,
			       ssize_t *idynamic, ssize_t *dim)
{
	(void)vars; // unused

	if (*peffects == effects)
		return;

	ssize_t delta = effects ? neffects : -neffects;

	*istatic += delta;
	*idynamic += delta;
	*dim += delta;
	*peffects = effects;
}

void design_set_send_effects(struct design *design, bool seffects)
{
	assert(design);

	design_set_effects(&design->seffects, seffects,
			   design_send_count(design), &design->send_vars,
			   &design->isstatic, &design->isdynamic,
			   &design->sdim);
}

void design_set_recv_effects(struct design *design, bool reffects)
{
	assert(design);

	design_set_effects(&design->reffects, reffects,
			   design_recv_count(design), &design->recv_vars,
			   &design->irstatic, &design->irdynamic,
			   &design->rdim);
}

static void design_add_var(struct design *design, const struct var_type *type,
			   void *params,
			   struct array *vars,
			   ssize_t *idynamic, ssize_t *ndynamic, ssize_t *dim)
{
	(void)idynamic; // unused
	struct design_var *var = array_add(vars, NULL);

	type->init(var, design, params);
	assert(var->dim >= 0);
	var->dyn_index = (*ndynamic);
	var->type = type;
	*ndynamic += var->dim;
	*dim += var->dim;
}

void design_add_send_var(struct design *design, const struct var_type *type,
			 void *params)
{
	assert(design);
	assert(type);
	assert(type->var_class == VAR_SEND_VAR);
	assert(type->init);

	design_add_var(design, type, params, &design->send_vars,
		       &design->isdynamic, &design->nsdynamic, &design->sdim);
}

void design_add_recv_var(struct design *design, const struct var_type *type,
			 void *params)
{
	assert(design);
	assert(type);
	assert(type->var_class == VAR_RECV_VAR);
	assert(type->init);

	design_add_var(design, type, params, &design->recv_vars,
		       &design->irdynamic, &design->nrdynamic, &design->rdim);
}

ssize_t design_send_traits_index(const struct design *design)
{
	assert(design);
	return design->isstatic;
}

ssize_t design_send_effects_index(const struct design *design)
{
	assert(design);

	if (design_send_effects(design))
		return design->iseffects;
	return -1;
}

ssize_t design_recv_effects_index(const struct design *design)
{
	assert(design);

	if (design_recv_effects(design))
		return design->ireffects;
	return -1;
}

static bool design_var_equals(const void *p1, const void *p2)
{
	const struct design_var *v1 = p1, *v2 = p2;
	return v1->type == v2->type;
}

static ssize_t design_var_dyn_index(const struct array *vars,
				    const struct var_type *type)
{
	assert(vars);
	assert(type);

	struct design_var *key = container_of(&type, struct design_var, type);

	const struct design_var *var =
	    array_find(vars, (predicate_fn) design_var_equals, key);

	if (var) {
		return var->dyn_index;
	}

	return -1;
}

ssize_t design_send_var_index(const struct design *design,
			      const struct var_type *type)
{
	assert(design);
	assert(type);

	ssize_t dyn_index = design_var_dyn_index(&design->send_vars, type);
	if (dyn_index >= 0) {
		ssize_t off = design_send_dyn_index(design);
		return off + dyn_index;
	} else {
		return -1;
	}
}

ssize_t design_recv_var_index(const struct design *design,
			      const struct var_type *type)
{
	assert(design);
	assert(type);

	ssize_t dyn_index = design_var_dyn_index(&design->recv_vars, type);
	if (dyn_index >= 0) {
		ssize_t off = design_recv_dyn_index(design);
		return off + dyn_index;
	} else {
		return -1;
	}

}
