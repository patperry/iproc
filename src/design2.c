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

#include "design2.h"


static void trait2_init_kvar(struct trait2 *v, const char *name,
			     const struct var *i, const struct var *j,
			     struct design2 *design, size_t index)
{
	assert(design_count(i->design) == design->count1);
	assert(design_count(j->design) == design->count2);
	assert(i->meta.type == VAR_TYPE_TRAIT);
	assert(j->meta.type == VAR_TYPE_TRAIT);

	size_t m = design->count1;
	size_t n = design->count2;
	size_t sizei = i->meta.size;
	size_t sizej = j->meta.size;
	size_t ranki = i->meta.rank;
	size_t rankj = j->meta.rank;
	size_t rank = ranki + rankj;
	size_t dims[VAR_RANK_MAX];
	memcpy(dims, i->meta.dims, ranki * sizeof(size_t));
	memcpy(dims + ranki, j->meta.dims, rankj * sizeof(size_t));

	var_meta_init(&v->var.meta, VAR_TYPE_TRAIT, name, dims, rank);
	v->var.design = design;
	v->var.index = index;

	v->dimi = sizei;
	v->xi = xmalloc(sizei * m * sizeof(double));
	if (sizei)
		lapack_dlacpy(LA_COPY_ALL, sizei, m,
			      design_trait_matrix(i->design) + i->index,
			      design_trait_dim(i->design), v->xi, sizei);

	v->dimj = sizej;
	v->xj = xmalloc(sizej * n * sizeof(double));
	if (sizej)
		lapack_dlacpy(LA_COPY_ALL, sizej, n,
			      design_trait_matrix(j->design) + j->index,
			      design_trait_dim(j->design), v->xj, sizej);

	v->xij = NULL;
}


static void trait2_deinit(struct trait2 *v)
{
	free(v->xij);
	free(v->xj);
	free(v->xi);
	var_meta_deinit(&v->var.meta);
}


static void tvar2_init(struct tvar2 *tv, const struct tvar2_type *type,
		       const char *name, struct design2 *design, size_t index,
		       va_list ap)
{
	size_t i, m = design2_count1(design);
	size_t n = design2_count2(design);

	type->init(&tv->var.meta, &tv->thunk, name, design, ap);
	tv->var.design = design;
	tv->var.index = index;
	tv->type = type;

	tv->deltaset = xmalloc(m * sizeof(struct deltaset));
	tv->tcur = xmalloc(m * sizeof(double));
	for (i = 0; i < m; i++) {
		deltaset_init(&tv->deltaset[i], n);
		tv->tcur[i] = -INFINITY;
	}
}


static void tvar2_clear(struct tvar2 *tv, size_t i)
{
	struct design2 *d = tv->var.design;
	size_t index = tv->var.index;
	size_t size = tv->var.meta.size;
	size_t dim = design2_tvar_dim(d);

	deltaset_clear(&tv->deltaset[i]);
	tv->tcur[i] = -INFINITY;

	double *ptr = d->tvar_x + d->tvar_ir[i] * dim + index;
	size_t iz, nz = d->tvar_ir[i+1] - d->tvar_ir[i];
	for (iz = 0; iz < nz; iz++) {
		memset(ptr, 0, size * sizeof(double));
		ptr += dim;
	}
}


static void tvar2_deinit(struct tvar2 *tv)
{
	struct design2 *d = tv->var.design;
	size_t i, m = design2_count1(d);

	for (i = 0; i < m; i++) {
		deltaset_deinit(&tv->deltaset[i]);
	}
	free(tv->deltaset);
	free(tv->tcur);

	tv->type->deinit(&tv->var.meta, tv->thunk, d);
}


void design2_init(struct design2 *d, struct history *h, size_t count1,
		size_t count2)
{
	assert(count2 == 0 || count1 <= SIZE_MAX / count2);  // ensure count1 * count2 < SIZE_MAX

	d->history = h;
	d->count1 = count1;
	d->count2 = count2;

	d->cohorts = xcalloc(count1, sizeof(*d->cohorts));
	d->cohort_reps = xcalloc(1, sizeof(*d->cohort_reps));
	d->ncohort = 1;

	d->trait_dim = 0;
	d->trait_x = NULL;
	d->traits = NULL;
	d->ntrait = 0;

	d->tvar_dim = 0;
	d->tvar_x = NULL;
	d->tvar_ir = xcalloc(count1 + 1, sizeof(size_t));
	d->tvar_jc = NULL;
	d->tvar_nnz = 0;
	d->tvar_nnz_max = 0;
	d->tvars = NULL;
	d->ntvar = 0;
	d->ntvar_max = 0;

	size_t i, m = count1;
	d->deltaset = xmalloc(m * sizeof(struct deltaset));
	d->history_version = xmalloc(m * sizeof(struct version_watch));
	d->tcur = xmalloc(m * sizeof(double));
	d->tnext = xmalloc(m * sizeof(double));

	for (i = 0; i < m; i++) {
		deltaset_init(&d->deltaset[i], count2);
		version_watch_init(&d->history_version[i], history_version(h));
		d->tcur[i] = -INFINITY;
		d->tnext[i] = INFINITY;
	}
}


static void design2_clear(struct design2 *d, size_t i)
{
	const struct history *h = design2_history(d);
	const struct version *v = history_version(h);

	size_t k, nk = d->ntvar;
	for (k = 0; k < nk; k++) {
		tvar2_clear(d->tvars[k], i);
	}

	deltaset_clear(&d->deltaset[i]);
	version_watch_set(&d->history_version[i], v);
	d->tcur[i] = -INFINITY;
	d->tnext[i] = INFINITY;
}


static void design2_tvars_update(struct design2 *d, size_t i)
{
	assert(i < design2_count1(d));
	const struct history *h = design2_history(d);
	const struct version *v = history_version(h);

	if (version_changed(v, &d->history_version[i])) {
		design2_clear(d, i);
	} else if (d->tcur[i] == history_time(h)) {
		return;
	}

	double t0 = d->tcur[i];
	d->tcur[i] = history_time(h);

	size_t k, nk = d->ntvar;
	struct tvar2 **tvars = d->tvars;
	double tnext = INFINITY;
	double t;

	for (k = 0; k < nk; k++) {
		if (tvars[k]->type->update) {
			t = tvars[k]->type->update(tvars[k], i, t0, h);
			tnext = MIN(tnext, t);
		}
	}

	d->tnext[i] = tnext;
}


const struct deltaset *design2_changes(const struct design2 *d, size_t i)
{
	assert(i < design2_count1(d));
	design2_tvars_update((struct design2 *)d, i);
	return &d->deltaset[i];
}


double design2_next_time(const struct design2 *d, size_t i)
{
	assert(i < design2_count1(d));
	design2_tvars_update((struct design2 *)d, i);
	return d->tnext[i];
}


static void traits_deinit(struct trait2 **traits, size_t len)
{
	struct trait2 *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = traits[i - 1];
		trait2_deinit(v);
	}
}


static void tvars_deinit(struct tvar2 **tvars, size_t len, struct design2 *d)
{
	struct tvar2 *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = tvars[i - 1];
		tvar2_deinit(v);
	}
}


void design2_deinit(struct design2 *d)
{
	struct version *hv = history_version(d->history);
	size_t i, m = design2_count1(d);

	for (i = 0; i < m; i++) {
		version_watch_deinit(&d->history_version[i], hv);
		deltaset_deinit(&d->deltaset[i]);
	}
	free(d->history_version);
	free(d->deltaset);

	tvars_deinit(d->tvars, d->ntvar, d);
	free2((void **)d->tvars, d->ntvar);
	free(d->tvar_jc);
	free(d->tvar_ir);
	free(d->tvar_x);

	traits_deinit(d->traits, d->ntrait);
	free2((void **)d->traits, d->ntrait);
	free(d->trait_x);

	free(d->cohort_reps);
	free(d->cohorts);
}


static void design2_traits_grow(struct design2 *d, size_t dntrait, size_t ddim)
{
	size_t ntrait = design2_trait_count(d);
	size_t ntrait1 = ntrait + dntrait;
	size_t dim = design2_trait_dim(d);
	size_t dim1 = dim + ddim;

	d->traits = xrealloc(d->traits, ntrait1 * sizeof(*d->traits));
	d->trait_dim = dim1;
}


static size_t cohort_strata_dim(struct design2 *d)
{
	size_t k, nk = d->ntrait;
	size_t dim = 0;

	for (k = 0; k < nk; k++) {
		const struct trait2 *v = d->traits[k];
		size_t vdim;

		if (v->xij) {
			vdim = v->var.meta.size * design2_count2(d);
		} else {
			vdim = v->dimi;
		}

		dim += vdim;
	}

	return dim;
}


static void recompute_cohorts(struct design2 *d)
{
	size_t dim = cohort_strata_dim(d);
	size_t k, nk = d->ntrait;
	size_t i, n = d->count1;

	struct strata s;
	double *x = xmalloc(dim * sizeof(double));

	strata_init(&s, dim);
	d->ncohort = 0;

	for (i = 0; i < n; i++) {
		size_t off = 0;

		for (k = 0; k < nk; k++) {
			const struct trait2 *v = d->traits[k];
			const double *vx;
			size_t vdim;

			if (v->xij) {
				vx = v->xij;
				vdim = v->var.meta.size * design2_count2(d);
			} else {
				vdim = v->dimi;
				vx = v->xi;
			}

			memcpy(x + off, vx + i * vdim, vdim * sizeof(double));
			off += vdim;
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
	size_t dim = d->trait_dim;
	size_t ic, nc = d->ncohort;
	size_t k, nk = d->ntrait;

	d->trait_x = xrealloc(d->trait_x, nc * n * dim * sizeof(double));
	memset(d->trait_x, 0, nc * n * dim * sizeof(double));

	// TODO: rewrite in terms of BLAS2

	for (ic = 0; ic < nc; ic++) {
		size_t i = d->cohort_reps[ic];
		for (j = 0; j < n; j++) {
			size_t off = 0;
			for (k = 0; k < nk; k++) {
				const struct trait2 *v = d->traits[k];
				double *dst = d->trait_x + ic * n * dim + j * dim + off;
				size_t size = v->var.meta.size;

				if (v->xij) {
					memcpy(dst, v->xij + i * n * size + j * size, size * sizeof(double));
				} else {
					const double *x = v->xj + j * v->dimj;
					const double *y = v->xi + i * v->dimi;

					if (v->dimj)
						blas_dger(v->dimj, v->dimi, 1.0, x, 1, y, 1, dst, v->dimj);
				}

				off += size;
			}
			assert(off == dim);
		}
	}
}


static struct trait2 *design2_trait_alloc_kron(struct design2 *d, size_t index,
					       const char *name, const struct var *i,
					       const struct var *j)
{
	assert(i->meta.rank + j->meta.rank <= VAR_RANK_MAX);
	struct trait2 *v = xmalloc(sizeof(*v));
	trait2_init_kvar(v, name, i, j, d, index);
	return v;
}


const struct var2 *design2_add_kron(struct design2 *d, const char *name,
				    const struct var *i, const struct var *j)
{
	assert(design_count(i->design) == d->count1);
	assert(design_count(j->design) == d->count2);
	assert(i->meta.type == VAR_TYPE_TRAIT); // other types not implemented
	assert(j->meta.type == VAR_TYPE_TRAIT); // ditto

	size_t dim0 = design2_trait_dim(d);
	size_t index = dim0;

	struct trait2 *v = design2_trait_alloc_kron(d, index, name, i, j);
	size_t size = v->var.meta.size;

	design2_traits_grow(d, 1, size);
	d->traits[d->ntrait++] = v;

	recompute_traits(d);

	return &v->var;
}


static void design2_grow_tvars(struct design2 *d, size_t delta)
{
	if (needs_grow(d->ntvar + delta, &d->ntvar_max)) {
		d->tvars = xrealloc(d->tvars, d->ntvar_max * sizeof(*d->tvars));
	}
}


const struct var2 *design2_add_tvar(struct design2 *d, const char *name, const struct tvar2_type *type, ...)
{
	assert(name);
	assert(type);

	struct tvar2 *tv = xmalloc(sizeof(*tv));
	struct var2 *v = &tv->var;
	va_list ap;

	va_start(ap, type);
	tvar2_init(tv, type, name, d, d->tvar_dim, ap);
	va_end(ap);

	size_t index = d->ntvar;
	design2_grow_tvars(d, 1);
	d->tvars[index] = tv;
	d->tvar_dim += v->meta.size;
	d->ntvar = index + 1;
	
	return v;
}


const struct var2 *design2_var(const struct design2 *d, const char *name)
{
	assert(name);

	size_t i, n;
	const struct var2 *v;

	n = d->ntrait;
	for (i = 0; i < n; i++) {
		v = &(d->traits[i]->var);
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


const double *design2_tvar(const struct design2 *d, const struct var2 *v, size_t i, size_t j)
{
	assert(v->design == d);
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	assert(v->meta.type == VAR_TYPE_TVAR);

	const double *x = design2_tvars(d, i, j);
	if (!x)
		return NULL;

	const size_t off = v->index;
	return x + off;
}


const double *design2_tvars(const struct design2 *d, size_t i, size_t j)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));

	design2_tvars_update((struct design2 *)d, i);

	size_t ir0 = d->tvar_ir[i];
	size_t ir1 = d->tvar_ir[i+1];
	size_t nz = ir1 - ir0;
	size_t *indx = d->tvar_jc + ir0;
	ptrdiff_t ix = find_index(j, indx, nz);

	if (ix < 0)
		return NULL;

	const double *x = d->tvar_x + (ir0 + ix) * d->tvar_dim;
	return x;
}


void design2_get_tvar_matrix(const struct design2 *d, size_t i, const double **x, const size_t **ind, size_t *nz)
{
	assert(i < design2_count1(d));
	
	design2_tvars_update((struct design2 *)d, i);

	size_t ir0 = d->tvar_ir[i];
	size_t ir1 = d->tvar_ir[i+1];

	*nz = ir1 - ir0;
	*ind = d->tvar_jc + ir0;
	*x = d->tvar_x + ir0 * d->tvar_dim;
}



double *design2_make_active(struct design2 *d, struct tvar2 *v, size_t i, size_t j)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));

	size_t dim = d->tvar_dim;
	size_t ir0 = d->tvar_ir[i];
	size_t ir1 = d->tvar_ir[i+1];
	size_t *indx = d->tvar_jc + ir0;
	size_t nz = ir1 - ir0;
	ptrdiff_t ix = find_index(j, indx, nz);

	if (ix < 0) {
		ix = ~ix;

		if (needs_grow(d->tvar_nnz + 1, &d->tvar_nnz_max)) {
			d->tvar_jc = xrealloc(d->tvar_jc, d->tvar_nnz_max * sizeof(size_t));
			d->tvar_x = xrealloc(d->tvar_x, d->tvar_nnz_max * dim * sizeof(double));
		}

		memmove(d->tvar_jc + (ir0 + ix + 1),
			d->tvar_jc + (ir0 + ix),
			(d->tvar_nnz - (ir0 + ix)) * sizeof(size_t));
		d->tvar_jc[ir0 + ix] = j;
		memmove(d->tvar_x + (ir0 + ix + 1) * dim,
			d->tvar_x + (ir0 + ix) * dim,
			(d->tvar_nnz - (ir0 + ix)) * dim * sizeof(double));
		memset(d->tvar_x + (ir0 + ix) * dim, 0, dim * sizeof(double));

		size_t i1, n1 = design2_count1(d);
		for (i1 = i + 1; i1 <= n1; i1++) {
			d->tvar_ir[i1] = d->tvar_ir[i1] + 1;
		}
		d->tvar_nnz = d->tvar_nnz + 1;

		assert(d->tvar_ir[n1] == d->tvar_nnz);
	}

	return d->tvar_x + (ir0 + ix) * dim + v->var.index;
}


void coefs2_init(struct coefs2 *c, const struct design2 *d)
{
	size_t dim0 = design2_trait_dim(d);
	size_t dim1 = design2_tvar_dim(d);
	c->traits = xcalloc(dim0, sizeof(double));
	c->tvars = xcalloc(dim1, sizeof(double));
}


void coefs2_deinit(struct coefs2 *c)
{
	free(c->tvars);
	free(c->traits);
}


void coefs2_set(struct coefs2 *dst, const struct coefs2 *src, const struct design2 *d)
{
	size_t dim0 = design2_trait_dim(d);
	size_t dim1 = design2_tvar_dim(d);

	if (src) {
		memcpy(dst->traits, src->traits, dim0 * sizeof(double));
		memcpy(dst->tvars, src->tvars, dim1 * sizeof(double));
	} else {
		memset(dst->traits, 0, dim0 * sizeof(double));
		memset(dst->tvars, 0, dim1 * sizeof(double));
	}
}

void coefs2_axpy(double alpha, const struct coefs2 *x, struct coefs2 *y, const struct design2 *d)
{
	size_t dim0 = design2_trait_dim(d);
	size_t dim1 = design2_tvar_dim(d);

	blas_daxpy(dim0, alpha, x->traits, 1, y->traits, 1);
	blas_daxpy(dim1, alpha, x->tvars, 1, y->tvars, 1);
}
