#include "port.h"
#include "design2.h"


static void prod2_init(struct tvar2 *tv, const char *name, struct history *h, va_list ap);
static void prod2_deinit(struct tvar2 *tv, struct history *h);
static void prod2_update_var(void *udata, struct design2 *d, const struct var2 *v, size_t i,
			     size_t j, const double *delta, const size_t *ind, size_t nz);

static struct design2_callbacks prod2_design_callbacks = {
	NULL, // update
	prod2_update_var,
	NULL // clear
};

static struct tvar2_type VAR2_PROD2_REP = {
	prod2_init,
	prod2_deinit
};

static const struct tvar2_type *VAR2_PROD2 = &VAR2_PROD2_REP;


const struct var2 *design2_add_prod(struct design2 *d, const char *name, const struct var2 *u, const struct var2 *v)
{
	assert(u->design == d);
	assert(v->design == d);
	assert(u->meta.size > 0);
	assert(v->meta.size > 0);

	const struct var2 *res = NULL;

	if (u->meta.type == VAR_TYPE_TRAIT && v->meta.type == VAR_TYPE_TRAIT) {
		assert(0 && "Not Implemented");
	} else {
		res = design2_add_tvar(d, name, VAR2_PROD2, d, u, v);
	}

	return res;
}


struct prod2_udata {
	struct design2 *design;
	const struct var2 *u;
	const struct var2 *v;
	double *delta;
	size_t *ind;
};

void prod2_init(struct tvar2 *tv, const char *name, struct history *h, va_list ap)
{
	struct design2 *d = va_arg(ap, struct design2*);
	const struct var2 *u = va_arg(ap, const struct var2*);
	const struct var2 *v = va_arg(ap, const struct var2*);

	size_t rank = u->meta.rank + v->meta.rank;
	size_t dims[VAR_RANK_MAX];
	memcpy(dims, u->meta.dims, u->meta.rank * sizeof(dims[0]));
	memcpy(dims + u->meta.rank, v->meta.dims, v->meta.rank * sizeof(dims[0]));
	var_meta_init(&tv->var.meta, name, VAR_TYPE_TVAR, dims, rank);
	assert(tv->var.meta.size == u->meta.size * v->meta.size);

	struct prod2_udata *udata = xmalloc(sizeof(*udata));
	udata->design = d;
	udata->u = u;
	udata->v = v;
	udata->delta = xmalloc(v->meta.size * sizeof(*udata->delta));
	udata->ind = xmalloc(v->meta.size * sizeof(*udata->ind));
	tv->udata = udata;

	design2_add_observer(d, tv, &prod2_design_callbacks);

	(void)h;
}


void prod2_deinit(struct tvar2 *tv, struct history *h)
{
	struct prod2_udata *udata = tv->udata;
	design2_remove_observer(udata->design, tv);
	free(udata->ind);
	free(udata->delta);
	free(udata);

	(void)h;
}


void prod2_update_var(void *udata, struct design2 *d, const struct var2 *v, size_t i,
		      size_t j, const double *delta, const size_t *ind, size_t nz)
{
	assert(ind || !nz);
	assert(nz <= v->meta.size);

	const struct tvar2 *tv = udata;
	const struct prod2_udata *udata0 = tv->udata;
	double *vdelta = udata0->delta;
	size_t *vind = udata0->ind;
	size_t vnz;
	size_t iz;

	if (!nz)
		return;

	if (v == udata0->u) {
		const double *dx = delta;
		const double *y;

		if (udata0->v->meta.type == VAR_TYPE_TRAIT) {
			y = design2_trait(d, udata0->v, i, j);
		} else {
			y = design2_tvar(d, udata0->v, i, j);
		}

		size_t iy, ny = udata0->v->meta.size;

		if (!y || !ny)
			return;

		vnz = nz * ny;
		memset(vdelta, 0, vnz * sizeof(*vdelta));
		blas_dger(ny, nz, 1.0, y, 1, dx, 1, vdelta, ny);

		size_t *dst = vind;

		for (iz = 0; iz < nz; iz++) {
			for (iy = 0; iy < ny; iy++) {
				*dst++ = ind[iz] * ny + iy;
			}
		}

	} else if (v == udata0->v) {
		const double *dy = delta;
		const double *x;

		if (udata0->u->meta.type == VAR_TYPE_TRAIT) {
			x = design2_trait(d, udata0->u, i, j);
		} else {
			x = design2_tvar(d, udata0->u, i, j);
		}

		size_t ix, nx = udata0->u->meta.size;
		size_t ny = udata0->v->meta.size;

		if (!x || !nx)
			return;

		vnz = nx * nz;
		memset(vdelta, 0, vnz * sizeof(*vdelta));
		blas_dger(nz, nx, 1.0, dy, 1, x, 1, vdelta, nz);

		size_t *dst = vind;

		for (ix = 0; ix < nx; ix++) {
			for (iz = 0; iz < nz; iz++) {
				*dst++ = ix * ny + ind[iz];
			}
		}
	} else {
		return;
	}

	design2_update(d, &tv->var, i, j, vdelta, vind, vnz);
}

