#ifndef DESIGN_H
#define DESIGN_H

#include "blas.h"
#include "sblas.h"
#include <stdarg.h>


struct design {
	struct frame *frame;
	size_t count;

	size_t ncohort;
	size_t *cohorts;
	size_t *cohort_reps;

	double *traits;
	struct var **trait_vars;
	size_t ntrait;

	size_t tvar_dim;
	struct tvar **tvars;
	size_t ntvar, ntvar_max;
	struct vpattern pat_buf;
	
	struct vpattern active;
	double *dx;	// dX[t]
	
	struct design_observer *observers;
	size_t nobs, nobs_max;
};

enum var_type {
	VAR_TYPE_TRAIT,
	VAR_TYPE_TVAR
};

struct var {
	const struct design *design;
	enum var_type type;
	const char *name;
	size_t dim;
	size_t index;
};

struct tvar {
	struct var var;
	const struct tvar_type *type;
	void *udata;
};

struct tvar_type {
	void (*init) (struct tvar *tv, const struct design *d, va_list ap);
	void (*deinit) (struct tvar * tv, const struct design *d);
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
};


void design_init(struct design *d, struct frame *f, size_t count);
void design_deinit(struct design *d);

/* observers */
void design_add_observer(struct design *d, void *udata,
			const struct design_callbacks *callbacks);
void design_remove_observer(struct design *d, void *udata);

/* properties */
static inline struct frame *design_frame(const struct design *d);
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
static inline const double *design_traits(const struct design *d);
const char *design_trait_name(const struct design *d, size_t k);
const struct var *design_add_trait(struct design *d, const char *name, const double *x);
void design_add_traits(struct design *d, const char * const *names, const double *x, size_t num);

void design_traits_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y);
void design_traits_tmul(double alpha, const struct design *d, const double *x, double beta, double *y);
void design_traits_axpy(double alpha, const struct design *d, size_t i, double *y);


/* tvars */
static inline size_t design_tvar_dim(const struct design *d);
static inline const double *design_tvars(const struct design *d, size_t i);
const char *design_tvar_name(const struct design *d, size_t k);
const struct var *design_add_tvar(struct design *d, const char *name, const struct tvar_type *type, ...);

static inline void design_tvars_get(const struct design *d, const double **dxp, const size_t **ip, size_t *nzp);
void design_tvar_get_lb(const struct design *d, size_t i, const double **dxp, const size_t **ip);
void design_tvar_get_ub(const struct design *d, size_t i, const double **dxp, const size_t **ip);

void design_tvars_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y);
void design_tvars_tmul(double alpha, const struct design *d, const double *x, double beta, double *y);
void design_tvars_axpy(double alpha, const struct design *d, size_t i, double *y);


/* coefs */
void coefs_init(struct coefs *c, const struct design *d);
void coefs_deinit(struct coefs *c);

void design_mul(double alpha, const struct design *d,
		const struct coefs *c, double beta, double *y);
void design_tmul(double alpha, const struct design *d, const double *x, double beta, struct coefs *c);
void design_axpy(double alpha, const struct design *d, size_t i, struct coefs *c);


/* internal functions (for use by tvar callbacks) */
void design_update(struct design *d, const struct var *v, size_t i, const double *delta,
		   const size_t *ind, size_t nz);


/* inline function definitions */
struct frame *design_frame(const struct design *d)
{
	return d->frame;
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
	return d->ntrait;
}


const double *design_traits(const struct design *d)
{
	return d->traits;
}


size_t design_tvar_dim(const struct design *d)
{
	return d->tvar_dim;
}


const double *design_tvars(const struct design *d, size_t i)
{
	assert(i < design_count(d));
	ptrdiff_t ix = vpattern_find(&d->active, i);

	if (ix < 0)
		return NULL;

	return d->dx + ix * d->tvar_dim;
}


void design_tvars_get(const struct design *d, const double **dxp, const size_t **ip, size_t *nzp)
{
	*dxp = d->dx;
	*ip = d->active.indx;
	*nzp = d->active.nz;
}


#endif /* DESIGN_H */
