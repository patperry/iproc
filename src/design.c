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

#include "vars.h"
#include "design.h"


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

static void design_frame_clear(void *udata, struct frame *f)
{
	struct design *d = udata;
	design_clear_range(d, 0, design_tvar_dim(d));
	
	size_t io, no = d->nobs;
	const struct design_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, d);
		}
	}
}

static struct frame_callbacks design_frame_callbacks = {
	NULL,
	NULL,
	design_frame_clear
};

void design_init(struct design *d, struct frame *f, size_t count)
{
	d->frame = f;
	d->count= count;

	d->cohorts = xcalloc(count, sizeof(*d->cohorts));
	d->cohort_reps = xcalloc(1, sizeof(*d->cohort_reps));
	d->ncohort = 1;	
	
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
	
	d->observers = NULL;
	d->nobs = 0;
	d->nobs_max = 0;
	
	frame_add_observer(f, d, &design_frame_callbacks);
}


static void tvars_deinit(struct tvar **tvars, size_t len, struct design *d)
{
	struct tvar *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = tvars[i - 1];
		if (v->type->deinit) {
			v->type->deinit(v, d);
		}
	}
}


void design_deinit(struct design *d)
{
	frame_remove_observer(d->frame, d);
	free(d->observers);
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
	size_t nmax = array_grow(d->nobs, d->nobs_max, delta, SIZE_MAX);
	if (nmax > d->nobs_max) {
		d->observers = xrealloc(d->observers, nmax * sizeof(d->observers[0]));
		d->nobs_max = nmax;
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


static void design_traits_grow(struct design *d, size_t delta)
{
	size_t count = design_count(d);		
	size_t dim = design_trait_dim(d);
	size_t dim1 = dim + delta;

	if (count) {
		double *traits = xmalloc(count * dim1 * sizeof(*traits));		
		lapack_dlacpy(LA_COPY_ALL, dim, count, d->traits, dim, traits, dim1);
		free(d->traits);
		d->traits = traits;
	}
	
	d->trait_vars = xrealloc(d->trait_vars, dim1 * sizeof(*d->trait_vars));
	d->ntrait = dim1;	
}


const char *design_trait_name(const struct design *d, size_t j)
{
	assert(j < design_trait_dim(d));
	return d->trait_vars[j]->name;
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

static struct var *design_trait_alloc(const struct design *d, const char *name,
				      size_t index, size_t dim)
{
	struct var *v = xmalloc(sizeof(*v));
	
	v->design = d;
	v->type = VAR_TYPE_TRAIT;
	v->name = xstrdup(name);
	v->dim = dim;
	v->index = index;
	
	return v;
}


const struct var *design_add_trait(struct design *d, const char *name, const double *x)
{
	size_t n = design_count(d);
	size_t ntrait = d->ntrait;
	size_t ntrait1 = ntrait + 1;
	struct var *v = design_trait_alloc(d, name, ntrait, 1);
	
	design_traits_grow(d, 1);
	d->trait_vars[ntrait] = v;
	blas_dcopy(n, x, 1, d->traits + ntrait, ntrait1);
	design_recompute_cohorts(d);
	return v;
}


void design_add_traits(struct design *d, const char * const *names, const double *x, size_t num)
{
	if (num == 0)
		return;

	size_t n = design_count(d);
	size_t ntrait = d->ntrait;
	size_t ntrait1 = ntrait;

	design_traits_grow(d, num);
	struct var *v;
	size_t i;
	
	for (i = 0; i < num; i++) {
		v = design_trait_alloc(d, names[i], ntrait1, 1);
		d->trait_vars[ntrait1] = v;
		ntrait1++;
	}
	assert(ntrait1 == ntrait + num);
	
	lapack_dlacpy(LA_COPY_ALL, num, n, x, num, d->traits + ntrait, ntrait1);
	design_recompute_cohorts(d);
}


const char *design_tvar_name(const struct design *d, size_t j)
{
	assert(j < design_tvar_dim(d));
	
	size_t i, n = d->ntvar;
	
	for (i = 0; i < n; i++) {
		const struct tvar *tv = d->tvars[i];
		const struct var *v = &tv->var;
		
		if (j < v->dim)
			return v->name;
		
		j -= v->dim;
	}
	
	assert(0);
	return NULL;
}


static void design_grow_tvars(struct design *d, size_t delta)
{
	size_t nmax = array_grow(d->ntvar, d->ntvar_max, delta, SIZE_MAX);
	if (nmax > d->ntvar_max) {
		d->tvars = xrealloc(d->tvars, nmax * sizeof(*d->tvars));
		d->ntvar_max = nmax;
	}
}


const struct var *design_add_tvar(struct design *d, const char *name, const struct tvar_type *type, ...)
{
	assert(name);
	assert(type);

	struct tvar *tv = xmalloc(sizeof(*tv));
	struct var *v = &tv->var;
	va_list ap;

	v->design = d;
	v->type = VAR_TYPE_TVAR;
	v->name = xstrdup(name);
	v->index = d->tvar_dim;
	v->dim = 0;
	tv->type = type;
	tv->udata = NULL;

	va_start(ap, type);
	type->init(tv, d, ap);
	va_end(ap);
	
	assert(tv->type == type);
	assert(v->index == d->tvar_dim);
	
	size_t index = d->ntvar;
	design_grow_tvars(d, 1);
	d->tvars[index] = tv;
	d->tvar_dim += v->dim;
	d->ind_buf = xrealloc(d->ind_buf, d->tvar_dim * sizeof(*d->ind_buf));
	d->ntvar = index + 1;
		
	return v;
}


void design_tvar_get_lb(const struct design *d, size_t i, const double **dxp, const size_t **ip)
{
	size_t ix = vpattern_lb(&d->active, i);
	*dxp = d->dx + ix * d->tvar_dim;
	*ip = d->active.indx + ix;
}


void design_tvar_get_ub(const struct design *d, size_t i, const double **dxp, const size_t **ip)
{
	size_t ix = vpattern_ub(&d->active, i);
	*dxp = d->dx + ix * d->tvar_dim;
	*ip = d->active.indx + ix;
}


const struct var *design_var(const struct design *d, const char *name)
{
	assert(name);

	size_t i, n;
	const struct var *v;

	n = d->ntrait;
	for (i = 0; i < n; i++) {
		v = d->trait_vars[i];
		if (strcmp(v->name, name) == 0) {
			return v;
		}
	}

	n = d->ntvar;
	for (i = 0; i < n; i++) {
		v = &(d->tvars[i]->var);
		if (strcmp(v->name, name) == 0) {
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
		assert(nz == v->dim);
		for (iz = 0; iz < v->dim; iz++) {
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
	assert(v->type == VAR_TYPE_TVAR);
	assert(i < design_count(d));
	
	size_t index = v->index;
	double *dx = design_dx(d, i) + index;
	
	if (ind) {
		sblas_daxpyi(nz, 1.0, delta, ind, dx);
	} else {
		assert(nz == v->dim);
		blas_daxpy(v->dim, 1.0, delta, 1, dx, 1);
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
	size_t dim0 = design_trait_dim(d);
	size_t dim1 = design_tvar_dim(d);
	size_t dim = dim0 + dim1;
	
	c->all = xmalloc(dim * sizeof(*c->all));
	c->traits = c->all;
	c->tvars = c->all + dim0;
	c->dim = dim;
}


void coefs_deinit(struct coefs *c)
{
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

