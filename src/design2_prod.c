#include "port.h"
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "xalloc.h"
#include "design2.h"



struct prod2_thunk {
	const struct var2 *u;
	const struct var2 *v;
};



static void prod2_init(struct var_meta *meta, void **thunk, const char *name, struct design2 *d, va_list ap)
{
	const struct var2 *u = va_arg(ap, const struct var2*);
	const struct var2 *v = va_arg(ap, const struct var2*);
	size_t rank = u->meta.rank + v->meta.rank;
	size_t dims[VAR_RANK_MAX];

	memcpy(dims, u->meta.dims, u->meta.rank * sizeof(size_t));
	memcpy(dims + u->meta.rank, v->meta.dims, v->meta.rank * sizeof(size_t));
	var_meta_init(meta, VAR_TYPE_TVAR, name, dims, rank);

	struct prod2_thunk *prod = xmalloc(sizeof(*prod));
	prod->u = u;
	prod->v = v;
	*thunk = prod;
}


static void prod2_deinit(struct var_meta *meta, void *thunk, struct design2 *d)
{
	struct prod2_thunk *prod = thunk;
	free(prod);
	var_meta_deinit(meta);
}


static double *get_x(const struct var2 *v, size_t i, size_t j)
{
	const struct design2 *d = v->design;
	size_t index = v->index;
	double *x;

	if (v->meta.type == VAR_TYPE_TRAIT) {
		size_t c = design2_cohort(d, i);
		size_t n = design2_count2(d);
		size_t dim = d->trait_dim;

		x = d->trait_x + (c * n + j) * dim + index;
	} else {
		size_t off = d->tvar_ir[i];
		const size_t *ind = d->tvar_jc + off;
		size_t nz = d->tvar_ir[i+1] - off;
		ptrdiff_t jz = find_index(j, ind, nz);
		size_t dim = d->tvar_dim;

		assert(v->meta.type == VAR_TYPE_TVAR);

		if (jz >= 0) {
			x = d->tvar_x + (off + jz) * dim + index;
		} else {
			x = NULL;
		}
	}

	return x;
}


static void design2_prod_update(struct design2 *d, struct tvar2 *prod_var,
				size_t i, size_t j, double t)
{
	struct prod2_thunk *prod = prod_var->thunk;
	const double *u = get_x(prod->u, i, j);
	const double *v = get_x(prod->v, i, j);
	size_t nu = prod->u->meta.size;
	size_t nv = prod->v->meta.size;
	double *x = design2_make_active(d, prod_var, i, j);
	size_t nx = prod_var->var.meta.size;

	memset(x, 0, nx * sizeof(double));

	if (u && v && nv) {
		blas_dger(nv, nu, 1.0, v, 1, u, 1, x, nv);
	}

	design2_update(d, prod_var, i, j, t);
}





static double prod2_update(struct tvar2 *prod_var, size_t i, double t0, const struct history *h)
{
	struct design2 *d = prod_var->var.design;
	struct prod2_thunk *prod = prod_var->thunk;
	struct delta *delta;

	if (prod->u->meta.type == VAR_TYPE_TVAR) {
		struct tvar2 *u = container_of(prod->u, struct tvar2, var);
		DELTASET_FOREACH(delta, &u->deltaset[i]) {
			double t = DELTA_TIME(delta);
			size_t j = DELTA_ITEM(delta);

			if (t < t0)
				break;

			design2_prod_update(d, prod_var, i, j, t);
		}
	}

	if (prod->v->meta.type == VAR_TYPE_TVAR) {
		struct tvar2 *v = container_of(prod->v, struct tvar2, var);
		DELTASET_FOREACH(delta, &v->deltaset[i]) {
			double t = DELTA_TIME(delta);
			size_t j = DELTA_ITEM(delta);

			if (t < t0)
				break;

			design2_prod_update(d, prod_var, i, j, t);
		}
	}

	return INFINITY; /* time of next update will be determined by u, v */
}


static struct tvar2_type VAR2_PROD2_REP = {
	prod2_init,
	prod2_deinit,
	prod2_update
};


const struct tvar2_type *VAR2_PROD2 = &VAR2_PROD2_REP;


static const struct var2 *add_trait_prod(struct design2 *d, const char *name,
					 const struct trait2 *u, const struct trait2 *v)
{
	const struct var2 *res = NULL;

	assert(0 && "Not Implemented");

	return res;
}


const struct var2 *design2_add_prod(struct design2 *d, const char *name,
				    const struct var2 *u, const struct var2 *v)
{
	assert(u->design == d);
	assert(v->design == d);

	const struct var2 *res = NULL;

	if (u->meta.type == VAR_TYPE_TRAIT && v->meta.type == VAR_TYPE_TRAIT) {
		res = add_trait_prod(d, name, container_of(u, struct trait2, var),
				     container_of(v, struct trait2, var));
	} else {
		res = design2_add_tvar(d, name, VAR2_PROD2, d, u, v);
	}

	return res;
}




