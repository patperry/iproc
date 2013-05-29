#ifndef DESIGN_H
#define DESIGN_H

#include "history.h"
#include "var.h"
#include "blas.h"
#include "deltaset.h"
#include "uintset.h"
#include <stdarg.h>


struct design_delta {
	struct uintset changed;
	int cleared;
};

struct design {
	struct history *history;
	size_t count;

	size_t ncohort;
	size_t *cohorts;
	size_t *cohort_reps;

	size_t trait_dim;
	double *traits;
	struct var **trait_vars;
	size_t ntrait;

	size_t tvar_dim;
	struct tvar **tvars;
	size_t ntvar, ntvar_max;
	size_t *ind_buf;
	
	struct uintset active;
	double *dx;	// dX[t]

	struct deltaset deltaset;

	struct design_observer *observers;
	size_t nobs, nobs_max;
};


struct design_callbacks {
	void (*update) (void *udata, struct design *d, size_t i,
			const double *delta, const size_t *ind, size_t nz);
	void (*update_var) (void *udata, struct design *d, const struct var *v, size_t i,
			    const double *delta, const size_t *ind, size_t nz);
	void (*clear) (void *udata, struct design *d);
};

struct design_observer {
	void *udata;
	struct design_callbacks callbacks;
};

struct coefs {
	double *all;
	double *traits; /* points into all */
	double *tvars;  /* points into all */
	size_t dim;
	int owner;
};


void design_init(struct design *d, struct history *h, size_t count);
void design_deinit(struct design *d);

/* observers */
void design_add_observer(struct design *d, void *udata,
			const struct design_callbacks *callbacks);
void design_remove_observer(struct design *d, void *udata);

/* properties */
static inline struct history *design_history(const struct design *d);
static inline size_t design_count(const struct design *d);
static inline size_t design_dim(const struct design *d);
const struct var *design_var(const struct design *d, const char *name);


/* cohorts */
static inline size_t design_cohort(const struct design *d, size_t i);
static inline size_t design_cohort_rep(const struct design *d, size_t c);
static inline size_t design_cohort_count(const struct design *d);
static inline void design_get_cohorts(const struct design *d,
				      const size_t **cohortsp,
				      const size_t **repsp, size_t *ncohortp);


/* traits */
static inline size_t design_trait_dim(const struct design *d);
static inline size_t design_trait_count(const struct design *d);
static inline const struct var * design_trait_var(const struct design *d, size_t k);
static inline const double *design_all_traits(const struct design *d);
static inline const double *design_traits(const struct design *d, size_t i);
static inline const double *design_trait(const struct design *d, const struct var *v, size_t i);

const struct var *design_add_trait(struct design *d, const char *name, const double *x, const size_t *dims, size_t rank);
void design_add_traits(struct design *d, const char * const *names, const double *x, size_t num);

void design_traits_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y);
void design_traits_tmul(double alpha, const struct design *d, const double *x, double beta, double *y);
void design_traits_axpy(double alpha, const struct design *d, size_t i, double *y);


/* tvars */
static inline size_t design_tvar_dim(const struct design *d);
static inline size_t design_tvar_count(const struct design *d);
static inline const struct var * design_tvar_var(const struct design *d, size_t k);
static inline const double *design_tvar(const struct design *d, const struct var *v, size_t i);
static inline const double *design_tvars(const struct design *d, size_t i);
const struct var *design_add_tvar(struct design *d, const char *name, const struct tvar_type *type, ...);

static inline void design_tvars_get_all(const struct design *d, const double **dxp, const size_t **ip, size_t *nzp);

void design_tvars_update(struct design *d);
void design_tvars_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y);
void design_tvars_tmul(double alpha, const struct design *d, const double *x, double beta, double *y);
void design_tvars_axpy(double alpha, const struct design *d, size_t i, double *y);


/* interactions */
const struct var *design_add_prod(struct design *d, const char *name, const struct var *u, const struct var *v);


/* coefs */
void coefs_init(struct coefs *c, const struct design *d);
void coefs_init_view(struct coefs *c, const struct design *d, const double *data);
void coefs_deinit(struct coefs *c);

void design_mul(double alpha, const struct design *d,
		const struct coefs *c, double beta, double *y);
void design_tmul(double alpha, const struct design *d, const double *x, double beta, struct coefs *c);
void design_axpy(double alpha, const struct design *d, size_t i, struct coefs *c);


/* internal functions (for use by tvar callbacks) */
void design_update(struct design *d, struct var *v, size_t i, const double *delta,
		   const size_t *ind, size_t nz);



/* inline function definitions */
struct history *design_history(const struct design *d)
{
	return d->history;
}


size_t design_count(const struct design *d)
{
	return d->count;
}


size_t design_dim(const struct design *d)
{
	return design_trait_dim(d) + design_tvar_dim(d);
}


size_t design_cohort(const struct design *d, size_t i)
{
	assert(i < design_count(d));
	return d->cohorts[i];
}


size_t design_cohort_rep(const struct design *d, size_t c)
{
	assert(c < design_cohort_count(d));
	return d->cohort_reps[c];
}


size_t design_cohort_count(const struct design *d)
{
	return d->ncohort;
}


void design_get_cohorts(const struct design *d, const size_t **cohortsp,
			const size_t **repsp, size_t *ncohortp)
{
	*cohortsp = d->cohorts;
	*repsp = d->cohort_reps;
	*ncohortp = d->ncohort;
}


size_t design_trait_dim(const struct design *d)
{
	return d->trait_dim;
}


size_t design_trait_count(const struct design *d)
{
	return d->ntrait;
}

const struct var * design_trait_var(const struct design *d, size_t k)
{
	assert(k < design_trait_count(d));
	return d->trait_vars[k];
}

const double *design_all_traits(const struct design *d)
{
	return d->traits;
}


const double *design_traits(const struct design *d, size_t i)
{
	assert(i < design_count(d));
	size_t dim = design_trait_dim(d);
	const double *x = design_all_traits(d);
	return x + i * dim;
}

const double *design_trait(const struct design *d, const struct var *v, size_t i)
{
	assert(v->design == d);
	assert(v->meta.type == VAR_TYPE_TRAIT);
	const double *x = design_traits(d, i);
	size_t off = v->index;
	return x + off;
}

size_t design_tvar_dim(const struct design *d)
{
	return d->tvar_dim;
}

size_t design_tvar_count(const struct design *d)
{
	return d->ntvar;
}

const struct var * design_tvar_var(const struct design *d, size_t k)
{
	assert(k <- design_tvar_count(d));
	return &d->tvars[k]->var;
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

	const double *dx = NULL;
	size_t iz;

	if (uintset_find(&d->active, i, &iz)) {
		dx = d->dx + iz * d->tvar_dim;
	}

	return dx;
}

void design_tvars_get_all(const struct design *d, const double **dxp, const size_t **ip, size_t *nzp)
{
	*dxp = d->dx;
	uintset_get_vals(&d->active, ip, nzp);
}


#endif /* DESIGN_H */
