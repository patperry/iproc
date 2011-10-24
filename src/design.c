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
	if (!strs)
		return NULL;

	char **res = xmalloc(len * sizeof(res[0]));
	size_t i;

	for (i = 0; i < len; i++) {
		res[i] = xstrdup(strs[i]);
	}
	return res;
}

void design_init(struct design *d, struct frame *f, size_t count)
{
	assert(d);

	d->frame = f;
	d->count= count;
	d->dim = 0;
	d->has_effects = 0;
	d->trait_off = 0;
	d->trait_dim = 0;
	d->traits = NULL;
	d->trait_names = NULL;
	d->dvar_off = 0;
	d->dvar_dim = 0;
	d->dvars = NULL;
	d->ndvar = 0;
	d->ndvar_max = 0;
}

static void free2(void **ptrs, size_t len)
{
	size_t i;

	if (!ptrs)
		return;

	for (i = len; i > 0; i--) {
		free(ptrs[i-1]);
	}
	free(ptrs);
}

static void free_dvars(struct design_var *dvars, size_t len, struct frame *f)
{
	struct design_var *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = &dvars[i - 1];
		frame_remove_observer(f, v);
		if (v->type->deinit) {
			v->type->deinit(v);
		}
	}
	free(dvars);
}

void design_deinit(struct design *d)
{
	assert(d);

	free_dvars(d->dvars, d->ndvar, d->frame);
	free2((void **)d->trait_names, d->trait_dim);
	free(d->traits);
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
	const struct vpattern *pat = &x->pattern;
	const size_t *indx = pat->indx;

	if (trans == BLAS_NOTRANS) {
		ptrdiff_t ix = vpattern_find(pat, off);
		size_t i = (ix < 0) ? ~ix : ix;

		for (; i < nz && indx[i] < off + dim; i++) {
                        blas_daxpy(count, alpha * dx[i],
				   a + (indx[i] - off) * lda, 1,
				   vector_to_ptr(y), 1);
                }
	} else {
		struct vector ysub = vector_slice(y, off, dim);

		sblas_dgemvi(trans, count, dim, alpha, a, lda,
			     dx, pat, 1.0, vector_to_ptr(&ysub));
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

void design_set_traits(struct design *d, const double *traits,
		       size_t dim, const char *const *names)
{
	free2((void **)d->trait_names, d->trait_dim);
	free(d->traits);

	if (dim >= d->trait_dim) {
		d->dim += (dim - d->trait_dim);
		d->dvar_off += (dim - d->trait_dim);
	} else {
		d->dim -= (d->trait_dim - dim);
		d->dvar_off -= (d->trait_dim - dim);
	}

	d->trait_dim = dim;
	d->traits = xmemdup(traits, d->count * dim * sizeof(traits[0]));
	d->trait_names = xstrdup2(names, dim);
}

static void design_grow_dvars(struct design *d)
{
	if (d->ndvar == d->ndvar_max) {
		struct frame *f = d->frame;
		size_t i;

		/* we need to update the frame observers since the call
		 * to realloc might relocate the dvars array
		 */
		for (i = d->ndvar; i > 0; i--) {
			frame_remove_observer(f, &d->dvars[i-1]);
		}

		size_t nmax = ARRAY_GROW(d->ndvar_max, SIZE_MAX);
		d->dvars = xrealloc(d->dvars, nmax * sizeof(d->dvars[0]));

		for (i = 0; i < d->ndvar; i++) {
			struct design_var *v = &d->dvars[i];
			frame_add_observer(f, v, &v->type->callbacks);
		}

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
	frame_add_observer(d->frame, v, &v->type->callbacks);

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

