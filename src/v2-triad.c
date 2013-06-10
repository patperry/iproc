#include "port.h"
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "xalloc.h"
#include "design2.h"



struct triad_thunk {
	const struct var2 *u;
	const struct var2 *v;
};



static void triad_init(struct var_meta *meta, void **thunk, const char *name, struct design2 *d, va_list ap)
{
	const struct var2 *u = va_arg(ap, const struct var2*);
	const struct var2 *v = va_arg(ap, const struct var2*);
	size_t rank = u->meta.rank + v->meta.rank;
	size_t dims[VAR_RANK_MAX];

	memcpy(dims, u->meta.dims, u->meta.rank * sizeof(size_t));
	memcpy(dims + u->meta.rank, v->meta.dims, v->meta.rank * sizeof(size_t));
	var_meta_init(meta, VAR_TYPE_TVAR, name, dims, rank);

	struct triad_thunk *triad = xmalloc(sizeof(*triad));
	triad->u = u;
	triad->v = v;
	*thunk = triad;
}


static void triad_deinit(struct var_meta *meta, void *thunk, struct design2 *d)
{
	struct triad_thunk *triad = thunk;
	free(triad);
	var_meta_deinit(meta);
}



static void design2_triad_recompute(struct design2 *d, struct tvar2 *triad_var,
				    size_t i)
{
	struct triad_thunk *triad = triad_var->thunk;
	
	assert(triad->u->meta.type == VAR_TYPE_TVAR); // not implemented otherwise
	assert(triad->v->meta.type == VAR_TYPE_TVAR); // ditto

	const size_t *ind = d->tvar_jc + d->tvar_ir[i];
	size_t iz, nz = d->tvar_ir[i+1] - d->tvar_ir[i];

	for (iz = 0; iz < nz; iz++) {
		size_t h = ind[iz];
		(void)h;
	}

	/*
	const double *u = get_x(triad->u, i, j);
	const double *v = get_x(triad->v, i, j);
	size_t nu = triad->u->meta.size;
	size_t nv = triad->v->meta.size;
	double *x = design2_make_active(d, triad_var, i, j);
	size_t nx = triad_var->var.meta.size;

	memset(x, 0, nx * sizeof(double));

	if (u && v && nv) {
		blas_dger(nv, nu, 1.0, v, 1, u, 1, x, nv);
	}

	design2_update(d, triad_var, i, j, t);
	 */
}


static double triad_update(struct tvar2 *triad_var, size_t i, double t0, const struct history *h)
{
	struct design2 *d = triad_var->var.design;
	struct triad_thunk *triad = triad_var->thunk;
	int update = 0;

	if (!design2_count2(d))
		return INFINITY;

	if (triad->u->meta.type == VAR_TYPE_TVAR) {
		struct tvar2 *u = container_of(triad->u, struct tvar2, var);
		if (deltaset_thead(&u->deltaset[i]) >= t0)
			update = 1;
	}

	if (!update && triad->v->meta.type == VAR_TYPE_TVAR) {
		struct tvar2 *v = container_of(triad->v, struct tvar2, var);
		const size_t *ind = d->tvar_jc + d->tvar_ir[i];
		size_t iz, nz = d->tvar_ir[i+1] - d->tvar_ir[i];

		for (iz = 0; iz < nz; iz++) {
			size_t j = ind[iz];
			if (deltaset_thead(&v->deltaset[j]) >= t0) {
				update = 1;
				break;
			}
		}
	}

	if (update)
		design2_triad_recompute(d, triad_var, i);

	return INFINITY; /* time of next update will be determined by u, v */
}


static struct tvar2_type VAR2_TRIAD_REP = {
	triad_init,
	triad_deinit,
	triad_update
};


const struct tvar2_type *VAR2_TRIAD = &VAR2_TRIAD_REP;
