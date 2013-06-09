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


static void trait_init(struct trait *v, const char *name, const size_t *dims,
		       size_t rank, struct design *design, size_t index)
{
	var_meta_init(&v->var.meta, VAR_TYPE_TRAIT, name, dims, rank);
	v->var.design = design;
	v->var.index = index;
}


static void trait_deinit(struct trait *v)
{
	var_meta_deinit(&v->var.meta);
}


static void tvar_init(struct tvar *tv, const struct tvar_type *type,
		      const char *name, struct design *design, size_t index, va_list ap)
{
	type->init(&tv->var.meta, &tv->thunk, name, design, ap);
	tv->var.design = design;
	tv->var.index = index;
	tv->type = type;
	deltaset_init(&tv->deltaset, design_count(design));
	tv->tcur = -INFINITY;
}


static void tvar_clear(struct tvar *tv)
{
	deltaset_clear(&tv->deltaset);
	tv->tcur = -INFINITY;

	struct design *d = tv->var.design;
	size_t index = tv->var.index;
	size_t size = tv->var.meta.size;

	size_t i, n = uintset_count(&d->active);
	size_t dim = design_tvar_dim(d);
	double *ptr = d->tvar_x + index;

	for (i = 0; i < n; i++) {
		memset(ptr, 0, size * sizeof(double));
		ptr += dim;
	}
}


static void tvar_deinit(struct tvar *tv)
{
	deltaset_deinit(&tv->deltaset);
	tv->type->deinit(&tv->var.meta, tv->thunk, tv->var.design);
}


void design_init(struct design *d, struct history *h, size_t count)
{
	d->history = h;
	d->count = count;

	d->cohorts = xcalloc(count, sizeof(*d->cohorts));
	d->cohort_reps = xcalloc(1, sizeof(*d->cohort_reps));
	d->ncohort = 1;	

	d->trait_dim = 0;
	d->trait_x = NULL;
	d->traits = NULL;
	d->ntrait = 0;

	d->tvar_dim = 0;
	d->tvar_x = NULL;
	uintset_init(&d->active);
	d->tvars = NULL;
	d->ntvar = 0;
	d->ntvar_max = 0;

	deltaset_init(&d->deltaset, count);
	version_watch_init(&d->history_version, history_version(h));
	d->tcur = -INFINITY;
	d->tnext = INFINITY;
}


static void design_clear(struct design *d)
{
	const struct history *h = design_history(d);
	const struct version *v = history_version(h);

	size_t i, n = d->ntvar;
	for (i = 0; i < n; i++) {
		tvar_clear(d->tvars[i]);
	}
	deltaset_clear(&d->deltaset);
	version_watch_set(&d->history_version, v);
	d->tcur = -INFINITY;
	d->tnext = INFINITY;
}


static void design_tvars_update(struct design *d)
{
	const struct history *h = design_history(d);
	const struct version *v = history_version(h);

	if (version_changed(v, &d->history_version)) {
		design_clear(d);
	} else if (d->tcur == history_time(h)) {
		return;
	}
	
	double t0 = d->tcur;
	d->tcur = history_time(h);

	size_t i, n = d->ntvar;
	struct tvar **tvars = d->tvars;
	double tnext = INFINITY;
	double t;

	for (i = 0; i < n; i++) {
		if (tvars[i]->type->update) {
			t = tvars[i]->type->update(tvars[i], t0, h);
			tnext = MIN(tnext, t);
		}
	}

	d->tnext = tnext;
}


const struct deltaset *design_changes(const struct design *d)
{
	design_tvars_update((struct design *)d);
	return &d->deltaset;
}


double design_next_time(const struct design *d)
{
	design_tvars_update((struct design *)d);
	return d->tnext;
}


static void traits_deinit(struct trait **traits, size_t len)
{
	struct trait *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = traits[i - 1];
		trait_deinit(v);
	}
}


static void tvars_deinit(struct tvar **tvars, size_t len)
{
	struct tvar *v;
	size_t i;

	for (i = len; i > 0; i--) {
		v = tvars[i - 1];
		tvar_deinit(v);
	}
}


void design_deinit(struct design *d)
{
	version_watch_deinit(&d->history_version,
			     history_version(d->history));
	deltaset_deinit(&d->deltaset);

	tvars_deinit(d->tvars, d->ntvar);
	free2((void **)d->tvars, d->ntvar);
	uintset_deinit(&d->active);
	free(d->tvar_x);

	traits_deinit(d->traits, d->ntrait);
	free2((void **)d->traits, d->ntrait);
	free(d->trait_x);

	free(d->cohort_reps);	
	free(d->cohorts);

}


static void design_traits_grow(struct design *d, size_t dntrait, size_t ddim)
{
	size_t count = design_count(d);

	size_t ntrait = design_trait_count(d);
	size_t ntrait1 = ntrait + dntrait;
	size_t dim = design_trait_dim(d);
	size_t dim1 = dim + ddim;

	if (count) {
		double *x = xmalloc(count * dim1 * sizeof(*x));
		lapack_dlacpy(LA_COPY_ALL, dim, count, d->trait_x, dim, x, dim1);
		free(d->trait_x);
		d->trait_x = x;
	}
	
	d->traits = xrealloc(d->traits, ntrait1 * sizeof(*d->traits));
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


static struct trait *design_trait_alloc(struct design *d, size_t index,
				      const char *name, const size_t *dims, size_t rank)
{
	assert(rank <= VAR_RANK_MAX);
	struct trait *v = xmalloc(sizeof(*v));
	trait_init(v, name, dims, rank, d, index);
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

	struct trait *v = design_trait_alloc(d, index, name, dims, rank);
	size_t size = v->var.meta.size;

	design_traits_grow(d, 1, size);
	d->traits[d->ntrait++] = v;

	size_t dim1 = dim0 + size;

	if (size) {
		lapack_dlacpy(LA_COPY_ALL, size, n, x, size, d->trait_x + dim0, dim1);
		design_recompute_cohorts(d);
	}
	
	return &v->var;
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
	struct trait *v;
	size_t i;
	
	for (i = 0; i < num; i++) {
		v = design_trait_alloc(d, index, names[i], NULL, 0);
		d->traits[ntrait1] = v;
		ntrait1++;
		index++;
	}
	d->ntrait = ntrait1;
	assert(d->trait_dim == dim1);

	lapack_dlacpy(LA_COPY_ALL, num, n, x, num, d->trait_x + dim0, dim1);
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
	assert(!d->tvar_x); // Not Implemented otherwise

	struct tvar *tv = xmalloc(sizeof(*tv));
	struct var *v = &tv->var;
	va_list ap;

	va_start(ap, type);
	tvar_init(tv, type, name, d, d->tvar_dim, ap);
	va_end(ap);

	size_t index = d->ntvar;
	design_grow_tvars(d, 1);
	d->tvars[index] = tv;
	d->tvar_dim += v->meta.size;
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
		v = &d->traits[i]->var;
		if (strcmp(v->meta.name, name) == 0) {
			return v;
		}
	}

	n = d->ntvar;
	for (i = 0; i < n; i++) {
		v = &d->tvars[i]->var;
		if (strcmp(v->meta.name, name) == 0) {
			return v;
		}
	}

	return NULL;
}


const double *design_tvar(const struct design *d, const struct var *v, size_t i)
{
	assert(v->design == d);
	assert(v->meta.type == VAR_TYPE_TVAR);

	const double *dx = design_tvars(d, i);
	if (!dx)
		return NULL;

	size_t off = v->index;
	return dx + off;
}


const double *design_tvars(const struct design *d, size_t i)
{
	assert(i < design_count(d));

	design_tvars_update((struct design *)d);

	const double *dx = NULL;
	size_t iz;

	if (uintset_find(&d->active, i, &iz)) {
		dx = d->tvar_x + iz * d->tvar_dim;
	}

	return dx;
}


void design_get_tvar_matrix(const struct design *d, const double **x, const size_t **ind, size_t *nz)
{
	design_tvars_update((struct design *)d);
	*x = d->tvar_x;
	uintset_get_vals(&d->active, ind, nz);
}


double *design_make_active(struct design *d, struct tvar *v, size_t i)
{
	size_t dim = d->tvar_dim;
	size_t iz;

	if (!uintset_find(&d->active, i, &iz)) {
		size_t nzmax0 = uintset_capacity(&d->active);
		size_t ntail = uintset_insert(&d->active, iz, i);
		size_t nzmax = uintset_capacity(&d->active);

		if (nzmax != nzmax0) {
			d->tvar_x = xrealloc(d->tvar_x, nzmax * dim * sizeof(*d->tvar_x));
		}

		memmove(d->tvar_x + (iz + 1) * dim,
			d->tvar_x + iz * dim,
			ntail * dim * sizeof(*d->tvar_x));
		memset(d->tvar_x + iz * dim, 0, dim * sizeof(*d->tvar_x));
	}

	return d->tvar_x + iz * dim + v->var.index;
}
