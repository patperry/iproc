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
#include "design2.h"


static void prod2_init(struct tvar2 *tv, struct design2 *d, va_list ap);
static void prod2_deinit(struct tvar2 *tv, struct design2 *d);
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



static void design2_clear_range(struct design2 *d, size_t joff, size_t len)
{
	assert(joff + len <= design2_tvar_dim(d));
	
	size_t i, n = d->nnz;
	size_t dim = design2_tvar_dim(d);
	double *ptr = d->dx + joff;
	
	for (i = 0; i < n; i++) {
		memset(ptr, 0, len * sizeof(*ptr));
		ptr += dim;
	}
}

static void design2_frame_clear(void *udata, struct frame *f)
{
	struct design2 *d = udata;
	design2_clear_range(d, 0, design2_tvar_dim(d));
	
	size_t io, no = d->nobs;
	const struct design2_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, d);
		}
	}
}




static void recompute_cohorts(struct design2 *d)
{

	size_t dim = d->kvar_dim;
	struct strata s;
	double *x = xmalloc(dim * sizeof(*x));

	strata_init(&s, dim);

	size_t i, n = d->count1;
	size_t k, nk = d->nkvar;

	d->ncohort = 0;

	for (i = 0; i < n; i++) {
		size_t off = 0;
		for (k = 0; k < nk; k++) {
			const struct kvar2 *kv = d->kvars[k];
			memcpy(x + off, kv->xi + i * kv->dimi, kv->dimi * sizeof(*x));
			off += kv->dimi;
		}
		assert(off == dim);

		size_t c = strata_add(&s, x);

		d->cohorts[i] = c;

		if (c == d->ncohort) {
			d->ncohort++;
			d->cohort_reps = xrealloc(d->cohort_reps, d->ncohort * sizeof(*d->cohort_reps));
			d->cohort_reps[c] = i;
		}

	}

	strata_deinit(&s);
	free(x);
}


static void recompute_traits(struct design2 *d)
{
	recompute_cohorts(d);

	size_t j, n = d->count2;
	size_t dim = d->kvar_dim;
	size_t ic, nc = d->ncohort;
	size_t k, nk = d->nkvar;

	d->traits = xrealloc(d->traits, nc * n * dim * sizeof(*d->traits));
	memset(d->traits, 0, nc * n * dim * sizeof(*d->traits));

	for (ic = 0; ic < nc; ic++) {
		size_t i = d->cohort_reps[ic];
		for (j = 0; j < n; j++) {
			size_t off = 0;
			for (k = 0; k < nk; k++) {
				const struct kvar2 *kv = d->kvars[k];
				const double *x = kv->xj + j * kv->dimj;
				const double *y = kv->xi + i * kv->dimi;
				double *a = d->traits + ic * n * dim + j * dim + off;

				if (kv->dimj)
					blas_dger(kv->dimj, kv->dimi, 1.0, x, 1, y, 1, a, kv->dimj);

				off += kv->dimi * kv->dimj;
			}
			assert(off == dim);
		}
	}

	d->ntrait = dim;
}


static struct frame_callbacks design2_frame_callbacks = {
	NULL,
	NULL,
	design2_frame_clear
};

void design2_init(struct design2 *d, struct frame *f, size_t count1, size_t count2)
{
	assert(count2 == 0 || count1 <= SIZE_MAX / count2);  // ensure count1 * count2 < SIZE_MAX

	d->frame = f;
	d->count1 = count1;
	d->count2 = count2;
	d->count = count1 * count2;

	d->ncohort = 1;
	d->cohorts = xcalloc(count1, sizeof(*d->cohorts));
	d->cohort_reps = xcalloc(1, sizeof(*d->cohort_reps));

	d->traits = NULL;
	d->ntrait = 0;

	d->kvar_dim = 0;
	d->kvars = NULL;
	d->nkvar = 0;
	d->nkvar_max = 0;

	d->tvar_dim = 0;
	d->tvars = NULL;
	d->ntvar = 0;
	d->ntvar_max = 0;
	d->ind_buf = NULL;
	
	d->ir = xcalloc(count1 + 1, sizeof(*d->ir));
	d->jc = NULL;
	d->dx = NULL;
	d->nnz = 0;
	d->nnz_max = 0;
	
	d->observers = NULL;
	d->nobs = 0;
	d->nobs_max = 0;
	
	frame_add_observer(f, d, &design2_frame_callbacks);
}


static void tvars_deinit(struct tvar2 **tvars, size_t len, struct design2 *d)
{
	struct tvar2 *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = tvars[i - 1];
		if (v->type->deinit) {
			v->type->deinit(v, d);
		}
	}
}

static void kvar_deinit(struct kvar2 *kvar)
{
	free(kvar->xi);
	free(kvar->xj);
}

static void kvars_deinit(struct kvar2 **kvars, size_t len)
{
	size_t i;
	struct kvar2 *v;

	for (i = len; i > 0; i--) {
		v = kvars[i - 1];
		kvar_deinit(v);
	}
}


void design2_deinit(struct design2 *d)
{
	frame_remove_observer(d->frame, d);
	free(d->observers);
	free(d->dx);
	free(d->jc);
	free(d->ir);
	free(d->ind_buf);
	tvars_deinit(d->tvars, d->ntvar, d);
	free2((void **)d->tvars, d->ntvar);
	kvars_deinit(d->kvars, d->nkvar);
	free2((void **)d->kvars, d->nkvar);
	free(d->traits);
}


static void design2_observers_grow(struct design2 *d, size_t delta)
{
	size_t nmax = array_grow(d->nobs, d->nobs_max, delta, SIZE_MAX);
	if (nmax > d->nobs_max) {
		d->observers = xrealloc(d->observers, nmax * sizeof(d->observers[0]));
		d->nobs_max = nmax;
	}
}


void design2_add_observer(struct design2 *d, void *udata,
			const struct design2_callbacks *callbacks)
{
	assert(d);
	assert(udata);
	assert(callbacks);
	
	design2_observers_grow(d, 1);
	struct design2_observer *obs = &d->observers[d->nobs++];
	obs->udata = udata;
	obs->callbacks = *callbacks;
}


void design2_remove_observer(struct design2 *d, void *udata)
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

/*
static void design2_traits_grow(struct design2 *d, size_t delta)
{
	size_t count = d->count;
	size_t dim = d->ntrait;
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


const char *design2_trait_name(const struct design2 *d, size_t j)
{
	assert(j < design2_trait_dim(d));
	return d->trait_vars[j]->name;
}


static struct var2 *design2_trait_alloc(const struct design2 *d, const char *name,
				      size_t index, size_t dim)
{
	struct var2 *v = xmalloc(sizeof(*v));
	
	v->design = d;
	v->type = VAR_TYPE_TRAIT;
	v->name = xstrdup(name);
	v->dim = dim;
	v->index = index;
	
	return v;
}


const struct var2 *design2_add_trait(struct design2 *d, const char *name, const double *x)
{
	size_t n = d->count;
	size_t ntrait = d->ntrait;
	size_t ntrait1 = ntrait + 1;
	struct var2 *v = design2_trait_alloc(d, name, ntrait, 1);
	
	design2_traits_grow(d, 1);
	d->trait_vars[ntrait] = v;
	blas_dcopy(n, x, 1, d->traits + ntrait, ntrait1);

	return v;
}


void design2_add_traits(struct design2 *d, const char * const *names, const double *x, size_t num)
{
	if (num == 0)
		return;
	
	size_t n = d->count;
	size_t ntrait = d->ntrait;
	size_t ntrait1 = ntrait;
	
	design2_traits_grow(d, num);
	struct var2 *v;
	size_t i;
	
	for (i = 0; i < num; i++) {
		v = design2_trait_alloc(d, names[i], ntrait1, 1);
		d->trait_vars[ntrait1] = v;
		ntrait1++;
	}
	assert(ntrait1 == ntrait + num);
	
	lapack_dlacpy(LA_COPY_ALL, num, n, x, num, d->traits + ntrait, ntrait1);
}
 */


static void design2_grow_kvars(struct design2 *d, size_t delta)
{
	size_t nmax = array_grow(d->nkvar, d->nkvar_max, delta, SIZE_MAX);
	if (nmax > d->nkvar_max) {
		d->kvars = xrealloc(d->kvars, nmax * sizeof(*d->kvars));
		d->nkvar_max = nmax;
	}
}

const struct var2 *design2_add_kron(struct design2 *d, const char *name,
				    const struct var *i, const struct var *j)
{
	assert(design_count(i->design) == d->count1);
	assert(design_count(j->design) == d->count2);
	assert(i->type == VAR_TYPE_TRAIT); // other types not implemented
	assert(j->type == VAR_TYPE_TRAIT); // ditto

	size_t m = d->count1;
	size_t n = d->count2;
	size_t dimi = i->dim;
	size_t dimj = j->dim;
	struct kvar2 *kv = xmalloc(sizeof(*kv));
	struct var2 *v = &kv->var;

	v->design = d;
	v->type = VAR_TYPE_TRAIT;
	v->name = xstrdup(name);
	v->index = d->kvar_dim;
	v->dim = dimi * dimj;

	kv->dimi = dimi;
	kv->xi = xmalloc(dimi * m * sizeof(*kv->xi));
	if (dimi)
		lapack_dlacpy(LA_COPY_ALL, dimi, m, design_all_traits(i->design) + i->index, dimi, kv->xi, dimi);

	kv->dimj = dimj;
	kv->xj = xmalloc(dimj * n * sizeof(*kv->xj));
	if (dimj)
		lapack_dlacpy(LA_COPY_ALL, dimj, n, design_all_traits(j->design) + j->index, dimj, kv->xj, dimj);

	size_t index = d->nkvar;
	design2_grow_kvars(d, 1);
	d->kvars[index] = kv;
	d->nkvar++;
	d->kvar_dim += v->dim;


	recompute_traits(d);

	return v;
}


const double *design2_all_traits(const struct design2 *d, size_t i)
{
	assert(i < design2_count1(d));
	size_t c = d->cohorts[i];
	size_t n = design2_count2(d);
	size_t dim = design2_trait_dim(d);
	const double *x = d->traits + c * n * dim;
	return x;
}


const double *design2_traits(const struct design2 *d, size_t i, size_t j)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	const double *x = design2_all_traits(d, i);
	size_t dim = design2_trait_dim(d);
	return x + j * dim;
}


const double *design2_trait(const struct design2 *d, const struct var2 *v, size_t i, size_t j)
{
	assert(v->design == d);
	assert(v->type == VAR_TYPE_TRAIT);
	const double *x = design2_traits(d, i, j);
	return x + v->index;
}


const char *design2_tvar_name(const struct design2 *d, size_t j)
{
	assert(j < design2_tvar_dim(d));
	
	size_t i, n = d->ntvar;
	
	for (i = 0; i < n; i++) {
		const struct tvar2 *tv = d->tvars[i];
		const struct var2 *v = &tv->var;
		
		if (j < v->dim)
			return v->name;
		
		j -= v->dim;
	}
	
	assert(0);
	return NULL;
}


static void design2_grow_tvars(struct design2 *d, size_t delta)
{
	size_t nmax = array_grow(d->ntvar, d->ntvar_max, delta, SIZE_MAX);
	if (nmax > d->ntvar_max) {
		d->tvars = xrealloc(d->tvars, nmax * sizeof(*d->tvars));
		d->ntvar_max = nmax;
	}
}


const struct var2 *design2_add_tvar(struct design2 *d, const char *name, const struct tvar2_type *type, ...)
{
	assert(name);
	assert(type);

	struct tvar2 *tv = xmalloc(sizeof(*tv));
	struct var2 *v = &tv->var;
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
	design2_grow_tvars(d, 1);
	d->tvars[index] = tv;
	d->tvar_dim += v->dim;
	d->ind_buf = xrealloc(d->ind_buf, d->tvar_dim * sizeof(*d->ind_buf));
	d->ntvar = index + 1;
	
	return v;
}


const struct var2 *design2_var(const struct design2 *d, const char *name)
{
	assert(name);

	size_t i, n;
	const struct var2 *v;

	//n = d->ntrait;
	//for (i = 0; i < n; i++) {
	//	v = d->trait_vars[i];
	//	if (strcmp(v->name, name) == 0) {
	//		return v;
	//	}
	//}

	n = d->nkvar;
	for (i = 0; i < n; i++) {
		v = &(d->kvars[i]->var);
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



void design2_traits_mul(double alpha, const struct design2 *d, size_t i,
			const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	size_t n = design2_count2(d);
	size_t p = design2_trait_dim(d);
	const double *a = design2_all_traits(d, i);
	size_t lda = p;	
	
	if (p) {
		blas_dgemv(BLAS_TRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	} else if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(n, beta, y, 1);
	}
}


void design2_traits_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));
	
	size_t n = design2_count2(d);
	size_t p = design2_trait_dim(d);
	const double *a = design2_all_traits(d, i);
	size_t lda = p;
	
	if (p) {
		blas_dgemv(BLAS_NOTRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	}

}


void design2_traits_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	
	size_t p = design2_trait_dim(d);
	const double *x = design2_traits(d, i, j);

	blas_daxpy(p, alpha, x, 1, y, 1);
}


void design2_tvars_mul(double alpha, const struct design2 *d, size_t i,
		       const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	const double *a;
	const size_t *j;
	size_t n = design2_count2(d);
	size_t dim = design2_tvar_dim(d);
	size_t nz;
	
	if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(n, beta, y, 1);
	}
	
	design2_tvars_get_all(d, i, &a, &j, &nz);
	for (; nz != 0; a += dim, j++, nz--) {
		y[*j] += alpha * blas_ddot(dim, a, 1, x, 1);
	}
}


void design2_tvars_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));
	
	size_t dim = design2_tvar_dim(d);
	const double *a;
	const size_t *j;
	size_t nz;
	
	if (beta == 0.0) {
		memset(y, 0, dim * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(dim, beta, y, 1);
	}
	
	design2_tvars_get_all(d, i, &a, &j, &nz);
	for (; nz != 0; a += dim, j++, nz--) {
		blas_daxpy(dim, alpha * x[*j], a, 1, y, 1);
	}
}


void design2_tvars_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	
	const double *dx = design2_tvars(d, i, j);
	
	if (dx) {
		size_t dim = design2_tvar_dim(d);
		blas_daxpy(dim, alpha, dx, 1, y, 1);
	}
}


static double *design2_dx(struct design2 *d, size_t i, size_t j)
{
	assert(i < design2_count1(d));

	size_t dim = d->tvar_dim;	
	size_t ir0 = d->ir[i];
	size_t ir1 = d->ir[i+1];
	size_t *indx = d->jc + ir0;
	size_t nz = ir1 - ir0;
	struct vpattern pat = vpattern_make(indx, nz);
	ptrdiff_t ix = vpattern_find(&pat, j);
	
	if (ix < 0) {
		ix = ~ix;
		
		if (d->nnz == d->nnz_max) {
			d->nnz_max = array_grow(d->nnz, d->nnz_max, 1, SIZE_MAX);
			d->jc = xrealloc(d->jc, d->nnz_max * sizeof(*d->jc));
			d->dx = xrealloc(d->dx, d->nnz_max * dim * sizeof(*d->dx));
		}
		
		assert(d->nnz < d->nnz_max);

		memmove(d->jc + (ir0 + ix + 1),
			d->jc + (ir0 + ix),
			(d->nnz - (ir0 + ix)) * sizeof(*d->jc));
		d->jc[ir0 + ix] = j;
		memmove(d->dx + (ir0 + ix + 1) * dim,
			d->dx + (ir0 + ix) * dim,
			(d->nnz - (ir0 + ix)) * dim * sizeof(*d->dx));
		memset(d->dx + (ir0 + ix) * dim, 0, dim * sizeof(*d->dx));
		
		size_t i1, n1 = design2_count1(d);
		for (i1 = i + 1; i1 <= n1; i1++) {
			d->ir[i1] = d->ir[i1] + 1;
		}
		d->nnz = d->nnz + 1;

		assert(d->ir[n1] == d->nnz);
	}
	
	return d->dx + (ir0 + ix) * dim;
}


static void design2_notify_update(struct design2 *d, const struct var2 *v,
				  size_t i, size_t j, const double *delta,
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
	const struct design2_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.update) {
			obs->callbacks.update(obs->udata, d, i, j, delta, d->ind_buf, nz);
		}
	}
}



void design2_update(struct design2 *d, const struct var2 *v, size_t i, size_t j, const double *delta,
		   const size_t *ind, size_t nz)
{
	assert(v->design == d);
	assert(v->type == VAR_TYPE_TVAR);
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	
	size_t index = v->index;
	double *dx = design2_dx(d, i, j) + index;
	
	if (ind) {
		sblas_daxpyi(nz, 1.0, delta, ind, dx);
	} else {
		assert(nz == v->dim);
		blas_daxpy(v->dim, 1.0, delta, 1, dx, 1);
	}
	
	size_t io, no = d->nobs;
	const struct design2_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.update_var) {
			obs->callbacks.update_var(obs->udata, d, v, i, j, delta, ind, nz);
		}
	}
	
	design2_notify_update(d, v, i, j, delta, ind, nz);
}


void coefs2_init(struct coefs2 *c, const struct design2 *d)
{
	size_t dim0 = design2_trait_dim(d);
	size_t dim1 = design2_tvar_dim(d);
	size_t dim = dim0 + dim1;
	
	c->all = xmalloc(dim * sizeof(*c->all));
	c->traits = c->all;
	c->tvars = c->all + dim0;
	c->dim = dim;
}


void coefs2_deinit(struct coefs2 *c)
{
	free(c->all);
}


void design2_mul(double alpha, const struct design2 *d, size_t i,
		const struct coefs2 *c, double beta, double *y)
{
	design2_traits_mul(alpha, d, i, c->traits, beta, y);
	design2_tvars_mul(alpha, d, i, c->tvars, 1.0, y);
}


void design2_tmul(double alpha, const struct design2 *d, size_t i,
		  const double *x, double beta, struct coefs2 *c)
{
	design2_traits_tmul(alpha, d, i, x, beta, c->traits);
	design2_tvars_tmul(alpha, d, i, x, beta, c->tvars);
}


void design2_axpy(double alpha, const struct design2 *d, size_t i, size_t j,
		  struct coefs2 *c)
{
	design2_traits_axpy(alpha, d, i, j, c->traits);
	design2_tvars_axpy(alpha, d, i, j, c->tvars);
}

const struct var2 *design2_add_prod(struct design2 *d, const char *name, const struct var2 *u, const struct var2 *v)
{
	assert(u->design == d);
	assert(v->design == d);
	assert(u->dim > 0);
	assert(v->dim > 0);

	const struct var2 *res = NULL;

	if (u->type == VAR_TYPE_TRAIT && v->type == VAR_TYPE_TRAIT) {
		assert(0 && "Not Implemented");
	} else {
		res = design2_add_tvar(d, name, VAR2_PROD2, u, v);
	}

	return res;
}


struct prod2_udata {
	const struct var2 *u;
	const struct var2 *v;
	double *delta;
	size_t *ind;
};

void prod2_init(struct tvar2 *tv, struct design2 *d, va_list ap)
{
	const struct var2 *u = va_arg(ap, const struct var2*);
	const struct var2 *v = va_arg(ap, const struct var2*);
	size_t dim = u->dim * v->dim;

	struct prod2_udata *udata = xmalloc(sizeof(*udata));
	udata->u = u;
	udata->v = v;
	udata->delta = xmalloc(dim * sizeof(*udata->delta));
	udata->ind = xmalloc(dim * sizeof(*udata->ind));

	tv->var.dim = dim;
	tv->udata = udata;

	design2_add_observer(d, tv, &prod2_design_callbacks);
}


void prod2_deinit(struct tvar2 *tv, struct design2 *d)
{
	design2_remove_observer(d, tv);

	struct prod2_udata *udata = tv->udata;
	free(udata->ind);
	free(udata->delta);
	free(udata);
}


void prod2_update_var(void *udata, struct design2 *d, const struct var2 *v, size_t i,
		      size_t j, const double *delta, const size_t *ind, size_t nz)
{
	assert(ind || !nz);
	assert(nz <= v->dim);

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

		if (udata0->v->type == VAR_TYPE_TRAIT) {
			y = design2_trait(d, udata0->v, i, j);
		} else {
			y = design2_tvar(d, udata0->v, i, j);
		}

		size_t iy, ny = udata0->v->dim;

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

		if (udata0->u->type == VAR_TYPE_TRAIT) {
			x = design2_trait(d, udata0->u, i, j);
		} else {
			x = design2_tvar(d, udata0->u, i, j);
		}

		size_t ix, nx = udata0->u->dim;
		size_t ny = udata0->v->dim;

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
