#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <search.h>
#include <stddef.h>
#include <stdlib.h>
#include "coreutil.h"
#include "sblas.h"
#include "xalloc.h"
#include "vars.h"
#include "design.h"

static char **xstrdup2(const char *const *strs, size_t len)
{
	char **res = xmalloc(len * sizeof(res[0]));
	size_t i;

	for (i = 0; i < len; i++) {
		res[i] = xstrdup(strs[i]);
	}
	return res;
}

void design_init(struct design *d, size_t count,
		 const double *intvls, size_t nintvl,
		 const double *traits, size_t trait_dim,
		 const char *const *trait_names)
{
	assert(d);
	assert(intvls || !nintvl);
#ifndef NDEBUG
	{
		size_t i;
		for (i = 1; i < nintvl; i++) {
			assert(intvls[i - 1] < intvls[i]);
		}
	}
#endif
	assert(traits || !trait_dim);

	d->count= count;
	d->intvls = xmemdup(intvls, nintvl * sizeof(intvls[0]));
	d->nintvl = nintvl;

	d->has_effects = 0;
	d->dim = trait_dim;
	d->trait_off = 0;
	d->trait_dim = trait_dim;
	d->traits = xmemdup(traits, count * trait_dim * sizeof(traits[0]));
	d->trait_names = xstrdup2(trait_names, trait_dim);
	d->dvar_off = trait_dim;
	d->dvar_dim = 0;
	d->dvars = NULL;
	d->ndvar = 0;
	d->ndvar_max = 0;
}

static void free2(void **ptrs, size_t len)
{
	size_t i;

	for (i = len; i > 0; i--) {
		free(ptrs[i-1]);
	}
	free(ptrs);
}

static void free_dvars(struct design_var *dvars, size_t len)
{
	struct design_var *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = &dvars[i - 1];
		if (v->type->deinit) {
			v->type->deinit(v);
		}
	}
	free(dvars);
}

void design_deinit(struct design *d)
{
	assert(d);

	free_dvars(d->dvars, d->ndvar);
	free2((void **)d->trait_names, d->trait_dim);
	free(d->traits);
	free(d->intvls);
}

static void
design_mul0_effects(double alpha,
		    enum blas_trans trans,
		    const struct design *d,
		    const struct vector *x, struct vector *y)
{
	if (!design_has_effects(d))
		return;

	size_t off = design_effects_index(d);
	size_t dim = design_count(d);

	if (trans == BLAS_NOTRANS) {
		struct vector xsub = vector_slice(x, off, dim);
		vector_axpy(alpha, &xsub, y);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		vector_axpy(alpha, x, &ysub);
	}
}

static void
design_muls0_effects(double alpha,
		     enum blas_trans trans,
		     const struct design *d,
		     const struct svector *x, struct vector *y)
{
	if (!design_has_effects(d))
		return;

	size_t off = design_effects_index(d);
	size_t dim = design_count(d);
	size_t end = off + dim;

	if (trans == BLAS_NOTRANS) {
		struct svector_iter itx;
		size_t i;
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
design_mul0_traits(double alpha,
		   enum blas_trans trans,
		   const struct design *d,
	 	   const struct vector *x, struct vector *y)
{
	if (!design_traits_dim(d))
		return;

	const double *traits = design_traits(d);
	size_t off = design_traits_index(d);
	size_t dim = design_traits_dim(d);
	size_t count = design_count(d);

	if (trans == BLAS_NOTRANS) {
		struct vector xsub = vector_slice(x, off, dim);
		blas_dgemv(trans, count, dim, alpha, traits, MAX(1, count),
			   vector_to_ptr(&xsub), 1, 1.0,
			   vector_to_ptr(y), 1);

	} else {
		struct vector ysub = vector_slice(y, off, dim);
		blas_dgemv(trans, count, dim, alpha, traits, MAX(1, count),
			   vector_to_ptr(x), 1, 1.0,
			   vector_to_ptr(&ysub), 1);
	}
}

static void
design_muls0_traits(double alpha,
		    enum blas_trans trans,
		    const struct design *d,
		    const struct svector *x, struct vector *y)
{
	if (!design_traits_dim(d))
		return;

	size_t off = design_traits_index(d);
	size_t dim = design_traits_dim(d);
	size_t count = design_count(d);
	const double *a = design_traits(d);
	size_t lda = MAX(1, count);
	size_t nz = svector_count(x);
	const double *dx = svector_data_ptr(x);
	const size_t *indx = (size_t *)svector_index_ptr(x);

	if (trans == BLAS_NOTRANS) {
		ptrdiff_t ix = sblas_find(nz, indx, off);
		size_t i = (ix < 0) ? ~ix : ix;

		for (; i < nz && indx[i] < off + dim; i++) {
                        blas_daxpy(count, alpha * dx[i],
				   a + (indx[i] - off) * lda, 1,
				   vector_to_ptr(y), 1);
                }
	} else {
		struct vector ysub = vector_slice(y, off, dim);

		sblas_dgemvi(trans, count, dim, nz, alpha, a, lda,
			     dx, indx, 1.0, vector_to_ptr(&ysub));
	}
}


void
design_mul0(double alpha, enum blas_trans trans, const struct design *d,
	    const struct vector *x, double beta, struct vector *y)
{
	assert(d);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || (size_t)vector_dim(x) == design_dim(d));
	assert(trans != BLAS_NOTRANS
	       || (size_t)vector_dim(y) == design_count(d));
	assert(trans == BLAS_NOTRANS
	       || (size_t)vector_dim(x) == design_count(d));
	assert(trans == BLAS_NOTRANS
	       || (size_t)vector_dim(y) == design_dim(d));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_mul0_effects(alpha, trans, d, x, y);
	design_mul0_traits(alpha, trans, d, x, y);
}


void
design_muls0(double alpha,
	     enum blas_trans trans,
	     const struct design *design,
	     const struct svector *x, double beta, struct vector *y)
{
	assert(design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || (size_t)svector_dim(x) == design_dim(design));
	assert(trans != BLAS_NOTRANS
	       || (size_t)vector_dim(y) == design_count(design));
	assert(trans == BLAS_NOTRANS
	       || (size_t)svector_dim(x) == design_count(design));
	assert(trans == BLAS_NOTRANS
	       || (size_t)vector_dim(y) == design_dim(design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	design_muls0_effects(alpha, trans, design, x, y);
	design_muls0_traits(alpha, trans, design, x, y);
}

void design_set_has_effects(struct design *d, int has_effects)
{
	assert(d);
	if (has_effects) {
		if (d->has_effects)
			return;

		d->dim += d->count;
		d->trait_off += d->count;
		d->dvar_off += d->count;
	} else {
		if (!d->has_effects )
			return;

		d->dim -= d->count;
		d->trait_off -= d->count;
		d->dvar_off -= d->count;
	}

	d->has_effects = has_effects;
}

static void design_grow_dvars(struct design *d)
{
	if (d->ndvar == d->ndvar_max) {
		size_t nmax = ARRAY_GROW(d->ndvar_max, SIZE_MAX);
		d->dvars = xrealloc(d->dvars, nmax * sizeof(d->dvars[0]));
		d->ndvar_max = nmax;
	}
}

void design_add_dvar(struct design *d, const struct var_type *type,
			 void *params)
{
	assert(d);
	assert(type);
	assert(type->init);

	design_grow_dvars(d);
	struct design_var *v = &d->dvars[d->ndvar];

	type->init(v, d, params);
	v->dyn_index = d->dvar_dim;
	v->type = type;

	d->ndvar++;
	d->dvar_dim += v->dim;
	d->dim += v->dim;
}

ptrdiff_t design_dvar_index(const struct design *d, const struct var_type *type)
{
	assert(d);
	assert(type);

	size_t i, n = d->ndvar;
	for (i = 0; i < n; i++) {
		struct design_var *v = &d->dvars[i];

		if (v->type == type) {
			return design_dvars_index(d) + v->dyn_index;
		}
	}

	return -1;
}

