#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <search.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "util.h"
#include "xalloc.h"

#include "vars.h"
#include "design2.h"


static void design2_clear_range(struct design2 *d, size_t joff, size_t len)
{
	assert(joff + len < design2_tvar_dim(d));
	
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
	
	d->traits = NULL;
	d->ldtraits = MAX(1, count1 * count2);
	d->trait_vars = NULL;
	d->ntrait = 0;
	d->ntrait_max = 0;

	d->tvar_dim = 0;
	d->tvars = NULL;
	d->ntvar = 0;
	d->ntvar_max = 0;
	d->pat_buf.indx = NULL;
	d->pat_buf.nzmax = 0;
	
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


void design2_deinit(struct design2 *d)
{
	frame_remove_observer(d->frame, d);
	free(d->observers);
	free(d->dx);
	free(d->jc);
	free(d->ir);
	free(d->pat_buf.indx);
	tvars_deinit(d->tvars, d->ntvar, d);
	free2((void **)d->tvars, d->ntvar);
	free2((void **)d->trait_vars, d->ntrait);
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


static void design2_traits_grow(struct design2 *d, size_t delta)
{
	size_t nmax = array_grow(d->ntrait, d->ntrait_max, delta, SIZE_MAX);
	if (nmax > d->ntrait_max) {
		size_t count = design2_count1(d) * design2_count2(d);
		d->traits = xrealloc(d->traits, nmax * count * sizeof(*d->traits));
		d->trait_vars = xrealloc(d->trait_vars, nmax * sizeof(*d->trait_vars));
		d->ntrait_max = nmax;
	}
}


const char *design2_trait_name(const struct design2 *d, size_t j)
{
	assert(j < design2_trait_dim(d));
	return d->trait_vars[j]->name;
}


const struct var2 *design2_add_trait(struct design2 *d, const char *name, const double *x)
{
	size_t ntrait = d->ntrait;
	struct var2 *v = xmalloc(sizeof(*v));
	
	v->design = d;
	v->type = VAR_TYPE_TRAIT;
	v->name = xstrdup(name);
	v->dim = 1;
	v->index = ntrait;
	
	design2_traits_grow(d, 1);
	d->trait_vars[ntrait] = v;
		
	memcpy(MATRIX_COL(&d->traits, d->ldtraits, ntrait), x, design2_count1(d) * design2_count2(d) * sizeof(*x));
	
	d->ntrait = ntrait + 1;

	return v;
}

void design2_add_traits(struct design2 *d, size_t ntrait, const char * const *names, const double *x)
{
	design2_traits_grow(d, ntrait);
	
	size_t ldx = d->ldtraits;
	size_t i;
	for (i = 0; i < ntrait; i++) {
		design2_add_trait(d, names[i], MATRIX_COL(x, ldx, i));
	}
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
	size_t nmax = array_grow(d->ntvar, d->ntrait_max, delta, SIZE_MAX);
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
	d->pat_buf.indx = xrealloc(d->pat_buf.indx, d->tvar_dim * sizeof(*d->pat_buf.indx));
	d->pat_buf.nzmax = d->tvar_dim;
	
	return v;
}


const struct var2 *design2_var(const struct design2 *d, const char *name)
{
	assert(name);

	size_t i, n;
	const struct var2 *v;

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



void design2_traits_mul(double alpha, const struct design2 *d, size_t i,
			const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	const double *a = design2_traits(d);
	size_t lda = d->ldtraits;
	size_t m = design2_count2(d);
	size_t n = design2_trait_dim(d);
	const double *asub = MATRIX_PTR(a, lda, i * m, 0);
	
	blas_dgemv(BLAS_NOTRANS, m, n, alpha, asub, lda, x, 1, beta, y, 1);
}


void design2_traits_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));
	
	const double *a = design2_traits(d);
	size_t lda = d->ldtraits;
	size_t m = design2_count2(d);
	size_t n = design2_trait_dim(d);
	const double *asub = MATRIX_PTR(a, lda, i * m, 0);
	
	blas_dgemv(BLAS_TRANS, m, n, alpha, asub, lda, x, 1, beta, y, 1);

}


void design2_traits_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	
	const double *a = design2_traits(d);
	size_t lda = d->ldtraits;
	size_t m = design2_count2(d);
	size_t dim = design2_trait_dim(d);

	blas_daxpy(dim, alpha, MATRIX_PTR(a, lda, i * m + j, 0), lda, y, 1);
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
	
	design2_tvars_get(d, i, &a, &j, &nz);
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
	
	design2_tvars_get(d, i, &a, &j, &nz);
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
				  const struct vpattern *pat)
{
	size_t index = v->index;
	size_t iz;
	
	if (!d->nobs)
		return;
	
	if (pat) {
		for (iz = 0; iz < pat->nz; iz++) {
			d->pat_buf.indx[iz] = pat->indx[iz] + index;
		}
		d->pat_buf.nz = pat->nz;
	} else {
		for (iz = 0; iz < v->dim; iz++) {
			d->pat_buf.indx[iz] = iz + index;
		}
		d->pat_buf.nz = v->dim;
	}
	
	size_t io, no = d->nobs;
	const struct design2_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.update) {
			obs->callbacks.update(obs->udata, d, i, j, delta, &d->pat_buf);
		}
	}
}



void design2_update(struct design2 *d, const struct var2 *v, size_t i, size_t j, const double *delta,
		   const struct vpattern *pat)
{
	assert(v->design == d);
	assert(v->type == VAR_TYPE_TVAR);
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	
	size_t index = v->index;
	double *dx = design2_dx(d, i, j) + index;
	
	if (pat) {
		sblas_daxpyi(1.0, delta, pat, dx);
	} else {
		blas_daxpy(v->dim, 1.0, delta, 1, dx, 1);
	}
	
	size_t io, no = d->nobs;
	const struct design2_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &d->observers[io];
		if (obs->callbacks.update_var) {
			obs->callbacks.update_var(obs->udata, d, v, i, j, delta, pat);
		}
	}
	
	design2_notify_update(d, v, i, j, delta, pat);
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
