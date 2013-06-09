#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <search.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "lapack.h"
#include "sblas.h"
#include "strata.h"
#include "util.h"
#include "xalloc.h"

#include "design.h"


struct prod_thunk {
	const struct var *u;
	const struct var *v;
};


static void prod_init(struct var_meta *meta, void **thunk, const char *name,
		      struct design *d, va_list ap)
{
	const struct var *u = va_arg(ap, const struct var*);
	const struct var *v = va_arg(ap, const struct var*);
	size_t rank = u->meta.rank + v->meta.rank;
	size_t dims[VAR_RANK_MAX];

	memcpy(dims, u->meta.dims, u->meta.rank * sizeof(size_t));
	memcpy(dims + u->meta.rank, v->meta.dims,
	       v->meta.rank * sizeof(size_t));
	var_meta_init(meta, VAR_TYPE_TVAR, name, dims, rank);

	struct prod_thunk *prod = xmalloc(sizeof(*prod));
	prod->u = u;
	prod->v = v;
	*thunk = prod;
}


static void prod_deinit(struct var_meta *meta, void *thunk, struct design *d)
{
	struct prod_thunk *prod = thunk;
	free(prod);
	var_meta_deinit(meta);
}


static double *get_x(const struct var *v, size_t i)
{
	const struct design *d = v->design;
	size_t index = v->index;
	size_t iz;
	double *x = NULL;

	if (v->meta.type == VAR_TYPE_TRAIT) {
		x = d->trait_x + i * d->trait_dim + index;
	} else {
		assert(v->meta.type == VAR_TYPE_TVAR);

		if (uintset_find(&d->active, i, &iz)) {
			x = d->tvar_x + iz * d->tvar_dim + index;
		}
	}

	return x;
}


static void design_prod_update(struct design *d, struct tvar *prod_var,
			       size_t i, double t)
{
	struct prod_thunk *prod = prod_var->thunk;
	const double *u = get_x(prod->u, i);
	const double *v = get_x(prod->v, i);
	size_t nu = prod->u->meta.size;
	size_t nv = prod->v->meta.size;
	double *x = design_make_active(d, prod_var, i);
	size_t nx = prod_var->var.meta.size;

	memset(x, 0, nx * sizeof(double));

	if (u && v && nv) {
		blas_dger(nv, nu, 1.0, v, 1, u, 1, x, nv);
	}

	design_update(d, prod_var, i, t);
}


static double prod_update(struct tvar *prod_var, double t0,
			  const struct history *h)
{
	struct design *d = prod_var->var.design;
	struct prod_thunk *prod = prod_var->thunk;
	struct delta *delta;

	if (prod->u->meta.type == VAR_TYPE_TVAR) {
		struct tvar *u = container_of(prod->u, struct tvar, var);
		DELTASET_FOREACH(delta, &u->deltaset) {
			double t = DELTA_TIME(delta);
			size_t i = DELTA_ITEM(delta);

			if (t < t0)
				break;

			design_prod_update(d, prod_var, i, t);
		}
	}

	if (prod->v->meta.type == VAR_TYPE_TVAR) {
		struct tvar *v = container_of(prod->v, struct tvar, var);
		DELTASET_FOREACH(delta, &v->deltaset) {
			double t = DELTA_TIME(delta);
			size_t i = DELTA_ITEM(delta);

			if (t < t0)
				break;

			design_prod_update(d, prod_var, i, t);
		}
	}

	return INFINITY; /* time of next update will be determined by u, v */
}


static struct tvar_type VAR_PROD_REP = {
	prod_init,
	prod_deinit,
	prod_update
};


const struct tvar_type *VAR_PROD = &VAR_PROD_REP;


const struct var *design_add_prod(struct design *d, const char *name,
				  const struct var *u, const struct var *v)
{
	assert(u->design == d);
	assert(v->design == d);

	const struct var *res = NULL;

	if (u->meta.type == VAR_TYPE_TRAIT && v->meta.type == VAR_TYPE_TRAIT) {
		size_t dims[VAR_RANK_MAX];
		size_t size = u->meta.size * v->meta.size;
		size_t rank = u->meta.rank + v->meta.rank;
		size_t i, n = design_count(d);
		double *x = xcalloc(n * size, sizeof(double));

		memcpy(dims, u->meta.dims, u->meta.rank * sizeof(size_t));
		memcpy(dims + u->meta.rank, v->meta.dims,
		       v->meta.rank * sizeof(size_t));

		for (i = 0; i < n; i++) {
			const double *xu = design_trait(d, u, i);
			const double *xv = design_trait(d, v, i);
			if (v->meta.size) {
				blas_dger(v->meta.size, u->meta.size, 1.0,
					  xv, 1, xu, 1, x + i * size,
					  v->meta.size);
			}
		}

		res = design_add_trait(d, name, x, dims, rank);
		free(x);
	} else {
		res = design_add_tvar(d, name, VAR_PROD, d, u, v);
	}

	return res;
}
