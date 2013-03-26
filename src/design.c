#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <search.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "lapack.h"
#include "strata.h"
#include "util.h"
#include "xalloc.h"

#include "design.h"


static void delta_init(struct design_delta *delta, size_t count)
{
	delta->cleared = 0;
	delta->ind = NULL;
	delta->dx = NULL;
	delta->nz = 0;
	delta->nzmax = 0;
}

static void delta_deinit(struct design_delta *delta)
{
	free(delta->dx);
	free(delta->ind);
}

void delta_clear(struct design_delta *delta)
{
	delta->nz = 0;
	delta->cleared = 0;
}

static void delta_set_cleared(struct design_delta *delta)
{
	delta->nz = 0;
	delta->cleared = 1;
}

static void delta_grow(struct design_delta *delta, size_t n, size_t dim)
{
	size_t nz1 = delta->nz + n;
	if (needs_grow(nz1, &delta->nzmax)) {
		delta->ind = xrealloc(delta->ind, delta->nzmax * sizeof(*delta->ind));
		delta->dx = xrealloc(delta->dx, delta->nzmax * dim * sizeof(*delta->dx));
	}
}

static size_t delta_search(struct design_delta *delta, size_t i, size_t dim)
{
	ptrdiff_t siz = find_index(i, delta->ind, delta->nz);
	size_t iz;

	/* not found */
	if (siz < 0) {
		iz = ~siz;
		size_t ntail = delta->nz - iz;

		delta_grow(delta, 1, dim);
		memmove(delta->ind + iz + 1, delta->ind + iz, ntail * sizeof(*delta->ind));
		memmove(delta->dx + (iz + 1) * dim, delta->dx + iz * dim, ntail * dim * sizeof(*delta->dx));

		delta->ind[iz] = i;
		memset(delta->dx + iz * dim, 0, dim * sizeof(*delta->dx));
		delta->nz++;
	} else {
		iz = siz;
	}

	return iz;
}

static void delta_update(struct design_delta *delta, const struct var *v, size_t i, const double *dx, const size_t *ind, size_t nz, size_t dim)
{
	size_t iz = delta_search(delta, i, dim);
	size_t index = v->index;
	double *dx_dst = delta->dx + iz * dim + index;

	if (ind) {
		sblas_daxpyi(nz, 1.0, dx, ind, dx_dst);
	} else {
		assert(nz == v->meta.size);
		blas_daxpy(v->meta.size, 1.0, dx, 1, dx_dst, 1);
	}
}



static void prod_init(struct tvar *tv, const char *name, struct history *h, va_list ap);
static void prod_deinit(struct tvar *tv, struct history *h);
static void prod_update_var(void *udata, struct design *d, const struct var *v, size_t i,
				  const double *delta, const size_t *ind, size_t nz);

static struct design_callbacks prod_design_callbacks = {
	NULL, // update
	prod_update_var,
	NULL // clear
};

static struct tvar_type VAR_PROD_REP = {
	prod_init,
	prod_deinit
};

static const struct tvar_type *VAR_PROD = &VAR_PROD_REP;



static void design_clear_range(struct design *d, size_t joff, size_t len)
{
	assert(joff + len <= design_tvar_dim(d));
	
	size_t i, n = d->active.nz;
	size_t dim = design_tvar_dim(d);
	double *ptr = d->dx + joff;
	
	for (i = 0; i < n; i++) {
		memset(ptr, 0, len * sizeof(*ptr));
		ptr += dim;
	}
}

static void design_history_clear(void *udata, struct history *h)
{
	struct design *d = udata;
	design_clear_range(d, 0, design_tvar_dim(d));

	size_t id, nd = d->ndelta;
	for (id = 0; id < nd; id++) {
		delta_set_cleared(d->deltas[id]);
	}

	size_t io, no = d->nobs;
	const struct design_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, d);
		}
	}

	(void)h;
}


static struct history_callbacks design_history_callbacks = {
       NULL,
       NULL,
       design_history_clear
};


void design_init(struct design *d, struct history *h, size_t count)
{
	d->history = h;
	d->count= count;

	d->cohorts = xcalloc(count, sizeof(*d->cohorts));
	d->cohort_reps = xcalloc(1, sizeof(*d->cohort_reps));
	d->ncohort = 1;	

	d->trait_dim = 0;
	d->traits = NULL;
	d->trait_vars = NULL;
	d->ntrait = 0;

	d->tvar_dim = 0;
	d->tvars = NULL;
	d->ntvar = 0;
	d->ntvar_max = 0;
	d->ind_buf = NULL;
	
	vpattern_init(&d->active);
	d->dx = NULL;

	d->deltas = NULL;
	d->ndelta = 0;
	d->ndelta_max = 0;

	d->observers = NULL;
	d->nobs = 0;
	d->nobs_max = 0;

	history_add_observer(h, d, &design_history_callbacks);
}


static void tvars_deinit(struct tvar **tvars, size_t len, struct design *d)
{
	struct history *h = design_history(d);
	struct tvar *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = tvars[i - 1];
		if (v->type->deinit) {
			v->type->deinit(v, h);
		}
	}
}


void design_deinit(struct design *d)
{
	history_remove_observer(d->history, d);
	free(d->observers);
	free(d->deltas);
	free(d->dx);
	vpattern_deinit(&d->active);
	free(d->ind_buf);
	tvars_deinit(d->tvars, d->ntvar, d);
	free2((void **)d->tvars, d->ntvar);
	free2((void **)d->trait_vars, d->ntrait);
	free(d->traits);
	free(d->cohort_reps);	
	free(d->cohorts);
}


static void design_observers_grow(struct design *d, size_t delta)
{
	if (needs_grow(d->nobs + delta, &d->nobs_max)) {
		d->observers = xrealloc(d->observers, d->nobs_max * sizeof(d->observers[0]));
	}
}


void design_add_observer(struct design *d, void *udata,
			const struct design_callbacks *callbacks)
{
	assert(d);
	assert(udata);
	assert(callbacks);
	
	design_observers_grow(d, 1);
	struct design_observer *obs = &d->observers[d->nobs++];
	obs->udata = udata;
	obs->callbacks = *callbacks;
}


void design_remove_observer(struct design *d, void *udata)
{
	assert(d);
	assert(udata);
	
	size_t i, n = d->nobs;
	for (i = n; i > 0; i--) {
		if (d->observers[i-1].udata == udata) {
			memmove(d->observers + i - 1, d->observers + i,
				(n - i) * sizeof(d->observers[0]));
			d->nobs = n - 1;
			break;
		}
	}
}


static void design_deltas_grow(struct design *d, size_t n)
{
	if (needs_grow(d->ndelta + n, &d->ndelta_max)) {
		d->deltas = xrealloc(d->deltas, d->ndelta_max * sizeof(d->deltas[0]));
	}
}


void design_delta_init(struct design *d, struct design_delta *delta)
{

	delta_init(delta, design_count(d));

	design_deltas_grow(d, 1);
	d->deltas[d->ndelta++] = delta;
}


void design_delta_deinit(struct design *d, struct design_delta *delta)
{
	size_t i, n = d->ndelta;
	for (i = n; i > 0; i--) {
		if (d->deltas[i-1] == delta) {
			memmove(d->deltas + i - 1, d->deltas + i,
				(n - i) * sizeof(d->deltas[0]));
			d->ndelta = n - 1;
			break;
		}
	}

	delta_deinit(delta);
}


void design_delta_clear(struct design *d, struct design_delta *delta)
{
	delta_clear(delta);
}


static void design_traits_grow(struct design *d, size_t dntrait, size_t ddim)
{
	size_t count = design_count(d);

	size_t ntrait = design_trait_count(d);
	size_t ntrait1 = ntrait + dntrait;
	size_t dim = design_trait_dim(d);
	size_t dim1 = dim + ddim;

	if (count) {
		double *traits = xmalloc(count * dim1 * sizeof(*traits));		
		lapack_dlacpy(LA_COPY_ALL, dim, count, d->traits, dim, traits, dim1);
		free(d->traits);
		d->traits = traits;
	}
	
	d->trait_vars = xrealloc(d->trait_vars, ntrait1 * sizeof(*d->trait_vars));
	d->trait_dim = dim1;
}


static void design_recompute_cohorts(struct design *d)
{
	size_t dim = design_trait_dim(d);
	size_t n = design_count(d);
	size_t i, c;
	
	struct strata strat;
	double *buf;
	
	strata_init(&strat, dim);
	buf = xmalloc(dim * sizeof(*buf));	

	d->ncohort = 0;
	for (i = 0; i < n; i++) {
		const double *x = design_traits(d, i);
		blas_dcopy(dim, x, 1, buf, 1);
		c = strata_add(&strat, buf);

		d->cohorts[i] = c;
		
		if (c == d->ncohort) {
			d->ncohort++;
			d->cohort_reps = xrealloc(d->cohort_reps, d->ncohort * sizeof(*d->cohort_reps));
			d->cohort_reps[c] = i;
		}
	}
	
	free(buf);
	strata_deinit(&strat);
}

static struct var *design_trait_alloc(struct design *d, const char *name,
				      size_t index, const size_t *dims, size_t rank)
{
	assert(rank <= VAR_RANK_MAX);
	struct var *v = xmalloc(sizeof(*v));
	
	var_meta_init(&v->meta, name, VAR_TYPE_TRAIT, dims, rank);
	v->design = d;
	v->index = index;
	
	return v;
}

const struct var *design_add_trait(struct design *d, const char *name,
				   const double *x, const size_t *dims,
				   size_t rank)
{
	assert(rank <= VAR_RANK_MAX);

	size_t n = design_count(d);
	size_t dim0 = design_trait_dim(d);
	size_t index = dim0;

	struct var *v = design_trait_alloc(d, name, index, dims, rank);
	size_t size = v->meta.size;


	design_traits_grow(d, 1, size);
	d->trait_vars[d->ntrait++] = v;

	size_t dim1 = dim0 + size;
	lapack_dlacpy(LA_COPY_ALL, size, n, x, size, d->traits + dim0, dim1);

	design_recompute_cohorts(d);
	return v;
}


void design_add_traits(struct design *d, const char * const *names, const double *x, size_t num)
{
	if (num == 0)
		return;

	size_t n = design_count(d);
	size_t dim0 = d->trait_dim;
	size_t dim1 = dim0 + num;
	size_t index = dim0;
	size_t ntrait = d->ntrait;
	size_t ntrait1 = ntrait;

	design_traits_grow(d, num, num);
	struct var *v;
	size_t i;
	
	for (i = 0; i < num; i++) {
		v = design_trait_alloc(d, names[i], index, NULL, 0);
		d->trait_vars[ntrait1] = v;
		ntrait1++;
		index++;
	}
	d->ntrait = ntrait1;
	assert(d->trait_dim == dim1);

	lapack_dlacpy(LA_COPY_ALL, num, n, x, num, d->traits + dim0, dim1);
	design_recompute_cohorts(d);
}


static void design_grow_tvars(struct design *d, size_t delta)
{
	if (needs_grow(d->ntvar + delta, &d->ntvar_max)) {
		d->tvars = xrealloc(d->tvars, d->ntvar_max * sizeof(*d->tvars));
	}
}


const struct var *design_add_tvar(struct design *d, const char *name, const struct tvar_type *type, ...)
{
	assert(name);
	assert(type);
	assert(!d->dx); // Not Implemented otherwise

	struct history *h = design_history(d);
	struct tvar *tv = xmalloc(sizeof(*tv));
	struct var *v = &tv->var;
	va_list ap;

	va_start(ap, type);
	type->init(tv, name, h, ap);
	va_end(ap);

	v->design = d;
	v->index = d->tvar_dim;
	tv->type = type;

	size_t index = d->ntvar;
	design_grow_tvars(d, 1);
	d->tvars[index] = tv;
	d->tvar_dim += v->meta.size;
	d->ind_buf = xrealloc(d->ind_buf, d->tvar_dim * sizeof(*d->ind_buf));
	d->ntvar = index + 1;
		
	return v;
}


const struct var *design_var(const struct design *d, const char *name)
{
	assert(name);

	size_t i, n;
	const struct var *v;

	n = d->ntrait;
	for (i = 0; i < n; i++) {
		v = d->trait_vars[i];
		if (strcmp(v->meta.name, name) == 0) {
			return v;
		}
	}

	n = d->ntvar;
	for (i = 0; i < n; i++) {
		v = &(d->tvars[i]->var);
		if (strcmp(v->meta.name, name) == 0) {
			return v;
		}
	}

	return NULL;
}


void design_traits_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y)
{
	size_t n = design_count(d);
	size_t p = design_trait_dim(d);
	const double *a = design_all_traits(d);
	size_t lda = p;

	if (p) {
		blas_dgemv(BLAS_TRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	} else if (beta == 0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1) {
		blas_dscal(n, beta, y, 1);
	}
}


void design_traits_tmul(double alpha, const struct design *d, const double *x, double beta, double *y)
{
	size_t n = design_count(d);
	size_t p = design_trait_dim(d);
	const double *a = design_all_traits(d);
	size_t lda = p;

	if (p) {
		blas_dgemv(BLAS_NOTRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	}
}


void design_traits_axpy(double alpha, const struct design *d, size_t i, double *y)
{
	assert(i < design_count(d));
	
	const double *x = design_traits(d, i);
	size_t dim = design_trait_dim(d);
	
	blas_daxpy(dim, alpha, x, 1, y, 1);
}


void design_tvars_mul(double alpha, const struct design *d,
		      const double *x, double beta, double *y)
{	
	size_t n = design_count(d);	
	size_t dim = design_tvar_dim(d);
	const double *a;
	const size_t *k;
	size_t nz;
	
	if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(n, beta, y, 1);
	}
	
	design_tvars_get_all(d, &a, &k, &nz);
	for (; nz != 0; a += dim, k++, nz--) {
		y[*k] += alpha * blas_ddot(dim, a, 1, x, 1);
	}
}


void design_tvars_tmul(double alpha, const struct design *d, const double *x, double beta, double *y)
{
	size_t dim = design_tvar_dim(d);
	const double *a;
	const size_t *k;
	size_t nz;
	
	if (beta == 0.0) {
		memset(y, 0, dim * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(dim, beta, y, 1);
	}

	design_tvars_get_all(d, &a, &k, &nz);
	for (; nz != 0; a += dim, k++, nz--) {
		blas_daxpy(dim, alpha * x[*k], a, 1, y, 1);
	}
}


void design_tvars_axpy(double alpha, const struct design *d, size_t i, double *y)
{
	assert(i < design_count(d));

	const double *dx = design_tvars(d, i);
	
	if (dx) {
		size_t dim = design_tvar_dim(d);
		blas_daxpy(dim, alpha, dx, 1, y, 1);
	}
}


static double *design_dx(struct design *d, size_t i)
{
	int ins;
	size_t dim = d->tvar_dim;
	size_t nzmax = d->active.nzmax;
	size_t ix = vpattern_search(&d->active, i, &ins);
	
	if (ins) {
		if (nzmax != d->active.nzmax) {
			nzmax = d->active.nzmax;
			d->dx = xrealloc(d->dx, nzmax * dim * sizeof(*d->dx));
		}
		
		memmove(d->dx + (ix + 1) * dim,
			d->dx + ix * dim,
			(d->active.nz - 1 - ix) * dim * sizeof(*d->dx));
		memset(d->dx + ix * dim, 0, dim * sizeof(*d->dx));
	}
	
	return d->dx + ix * dim;
}


static void design_notify_update(struct design *d, const struct var *v,
				 size_t i, const double *delta,
				 const size_t *ind, size_t nz)
{
	size_t index = v->index;
	size_t iz;
	
	if (!d->nobs)
		return;
	
	if (ind) {
		for (iz = 0; iz < nz; iz++) {
			d->ind_buf[iz] = ind[iz] + index;
		}
	} else {
		assert(nz == v->meta.size);
		for (iz = 0; iz < v->meta.size; iz++) {
			d->ind_buf[iz] = iz + index;
		}
	}
	
	size_t io, no = d->nobs;
	const struct design_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.update) {
			obs->callbacks.update(obs->udata, d, i, delta, d->ind_buf, nz);
		}
	}
}


void design_update(struct design *d, const struct var *v, size_t i, const double *delta,
		   const size_t *ind, size_t nz)
{
	assert(v->design == d);
	assert(v->meta.type == VAR_TYPE_TVAR);
	assert(i < design_count(d));
	
	size_t index = v->index;
	double *dx = design_dx(d, i) + index;
	
	if (ind) {
		sblas_daxpyi(nz, 1.0, delta, ind, dx);
	} else {
		assert(nz == v->meta.size);
		blas_daxpy(v->meta.size, 1.0, delta, 1, dx, 1);
	}

	size_t id, nd = d->ndelta;
	size_t dim = design_dim(d);
	for (id = 0; id < nd; id++) {
		delta_update(d->deltas[id], v, i, delta, ind, nz, dim);
	}

	size_t io, no = d->nobs;
	const struct design_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.update_var) {
			obs->callbacks.update_var(obs->udata, d, v, i, delta, ind, nz);
		}
	}

	design_notify_update(d, v, i, delta, ind, nz);
}


void coefs_init(struct coefs *c, const struct design *d)
{
	size_t dim = design_dim(d);
	double *data = xmalloc(dim * sizeof(data[0]));
	coefs_init_view(c, d, data);
	c->owner = 1;
}

void coefs_init_view(struct coefs *c, const struct design *d, const double *data)
{
	size_t dim0 = design_trait_dim(d);
	size_t dim1 = design_tvar_dim(d);
	size_t dim = dim0 + dim1;

	c->all = (double *)data;
	c->traits = c->all;
	c->tvars = c->all + dim0;
	c->dim = dim;
	c->owner = 0;
}


void coefs_deinit(struct coefs *c)
{
	if (c->owner)
		free(c->all);
}


void design_mul(double alpha, const struct design *d,
		const struct coefs *c, double beta, double *y)
{
	design_traits_mul(alpha, d, c->traits, beta, y);
	design_tvars_mul(alpha, d, c->tvars, 1.0, y);
}


void design_tmul(double alpha, const struct design *d, const double *x,
		 double beta, struct coefs *c)
{
	design_traits_tmul(alpha, d, x, beta, c->traits);
	design_tvars_tmul(alpha, d, x, beta, c->tvars);
}


void design_axpy(double alpha, const struct design *d, size_t i, struct coefs *c)
{
	design_traits_axpy(alpha, d, i, c->traits);
	design_tvars_axpy(alpha, d, i, c->tvars);
}

const struct var *design_add_prod(struct design *d, const char *name, const struct var *u, const struct var *v)
{
	assert(u->design == d);
	assert(v->design == d);
	assert(u->meta.size > 0);
	assert(v->meta.size > 0);

	const struct var *res = NULL;
	size_t size = u->meta.size * v->meta.size;
	size_t rank = u->meta.rank + v->meta.rank;
	size_t dims[VAR_RANK_MAX];
	size_t i, n = design_count(d);

	memcpy(dims, u->meta.dims, u->meta.rank * sizeof(dims[0]));
	memcpy(dims + u->meta.rank, v->meta.dims, v->meta.rank * sizeof(dims[0]));

	if (u->meta.type == VAR_TYPE_TRAIT && v->meta.type == VAR_TYPE_TRAIT) {
		double *x = xcalloc(n * size, sizeof(double));

		for (i = 0; i < n; i++) {
			const double *xu = design_trait(d, u, i);
			const double *xv = design_trait(d, v, i);
			blas_dger(v->meta.size, u->meta.size, 1.0, xv, 1, xu, 1,
				  x + i * size, v->meta.size);
		}

		res = design_add_trait(d, name, x, dims, rank);
		free(x);
	} else {
		res = design_add_tvar(d, name, VAR_PROD, d, u, v);
	}

	return res;
}



struct prod_udata {
	struct design *design;
	const struct var *u;
	const struct var *v;
	double *delta;
	size_t *ind;
};

void prod_init(struct tvar *tv, const char *name, struct history *h, va_list ap)
{
	struct design *d = va_arg(ap, struct design*);
	const struct var *u = va_arg(ap, const struct var*);
	const struct var *v = va_arg(ap, const struct var*);
	size_t size = u->meta.size * v->meta.size;
	size_t rank = u->meta.rank + v->meta.rank;
	size_t dims[VAR_RANK_MAX];

	memcpy(dims, u->meta.dims, u->meta.rank * sizeof(dims[0]));
	memcpy(dims + u->meta.rank, v->meta.dims, v->meta.rank * sizeof(dims[0]));
	var_meta_init(&tv->var.meta, name, VAR_TYPE_TVAR, dims, rank);

	struct prod_udata *udata = xmalloc(sizeof(*udata));
	udata->design = d;
	udata->u = u;
	udata->v = v;
	udata->delta = xmalloc(size * sizeof(*udata->delta));
	udata->ind = xmalloc(size * sizeof(*udata->ind));
	tv->udata = udata;

	design_add_observer(d, tv, &prod_design_callbacks);
}


void prod_deinit(struct tvar *tv, struct history *h)
{
	struct prod_udata *udata = tv->udata;
	design_remove_observer(udata->design, tv);
	free(udata->ind);
	free(udata->delta);
	free(udata);
}


void prod_update_var(void *udata, struct design *d, const struct var *v, size_t i,
		     const double *delta, const size_t *ind, size_t nz)
{
	assert(ind || !nz);
	assert(nz <= v->meta.size);

	const struct tvar *tv = udata;
	const struct prod_udata *udata0 = tv->udata;
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
			y = design_trait(d, udata0->v, i);
		} else {
			y = design_tvar(d, udata0->v, i);
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
			x = design_trait(d, udata0->u, i);
		} else {
			x = design_tvar(d, udata0->u, i);
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

	design_update(d, &tv->var, i, vdelta, vind, vnz);
}