#ifndef DESIGN_H
#define DESIGN_H

#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include "blas.h"
#include "deltaset.h"
#include "uintset.h"
#include "history.h"
#include "var.h"


struct var {
	struct var_meta meta;
	struct design *design;
	size_t index;
};

struct trait {
	struct var var;
};

struct tvar {
	struct var var;
	void *thunk;
	const struct tvar_type *type;
	struct deltaset deltaset;
	double tcur;
};


struct tvar_type {
	void (*init) (struct var_meta *meta, void **thunk, const char *name, struct design *d, va_list ap);
	void (*deinit) (struct var_meta *meta, void *thunk, struct design *d);
	void (*update) (struct tvar *tv, double t0, const struct history *h);
};

extern const struct tvar_type *VAR_IRECVTOT;
extern const struct tvar_type *VAR_ISENDTOT;
extern const struct tvar_type *VAR_NRECVTOT;
extern const struct tvar_type *VAR_NSENDTOT;


struct design {
	struct history *history;
	size_t count;

	size_t ncohort;
	size_t *cohorts;
	size_t *cohort_reps;

	size_t trait_dim;
	double *trait_x;
	struct trait **traits;
	size_t ntrait;

	size_t tvar_dim;
	double *tvar_x;
	struct uintset active;
	struct tvar **tvars;
	size_t ntvar, ntvar_max;

	struct deltaset deltaset;
	double tcur;
};


struct coefs {
	double *all;
	double *traits; /* points into all */
	double *tvars;  /* points into all */
	size_t dim;
	int owner;
};


/* create/destroy */
void design_init(struct design *d, struct history *h, size_t count);
void design_deinit(struct design *d);


/* properties */
static inline struct history *design_history(const struct design *d) { return d->history; }
static inline size_t design_count(const struct design *d) { return d->count; }
static inline size_t design_dim(const struct design *d) { return d->trait_dim + d->tvar_dim; }

const struct var *design_var(const struct design *d, const char *name);


/* cohorts */
static inline size_t design_cohort_count(const struct design *d) { return d->ncohort; }

static inline size_t design_cohort(const struct design *d, size_t i)
{
	assert(i < design_count(d));
	return d->cohorts[i];
}

static inline size_t design_cohort_rep(const struct design *d, size_t c)
{
	assert(c < design_cohort_count(d));
	return d->cohort_reps[c];
}

static inline void design_get_cohorts(const struct design *d, const size_t **cohorts,
				      const size_t **reps, size_t *ncohort)
{
	*cohorts = d->cohorts;
	*reps = d->cohort_reps;
	*ncohort = d->ncohort;
}


/* traits */
static inline size_t design_trait_count(const struct design *d) { return d->ntrait; }
static inline size_t design_trait_dim(const struct design *d) { return d->trait_dim; }

static inline const struct var *design_trait_item(const struct design *d, size_t k)
{
	assert(k < design_trait_count(d));
	return &d->traits[k]->var;
}

const struct var *design_add_trait(struct design *d, const char *name, const double *x,
				   const size_t *dims, size_t rank);
void design_add_traits(struct design *d, const char * const *names, const double *x, size_t num);


static inline const double *design_trait_matrix(const struct design *d)
{
	return d->trait_x;
}

static inline const double *design_traits(const struct design *d, size_t i)
{
	assert(i < design_count(d));
	size_t dim = design_trait_dim(d);
	const double *x = design_trait_matrix(d);
	return x + i * dim;
}

static inline const double *design_trait(const struct design *d, const struct var *v, size_t i)
{
	assert(v->design == d);
	assert(v->meta.type == VAR_TYPE_TRAIT);
	const double *x = design_traits(d, i);
	size_t off = v->index;
	return x + off;
}


/* tvars */
static inline size_t design_tvar_count(const struct design *d) { return d->ntvar; }
static inline size_t design_tvar_dim(const struct design *d) { return d->tvar_dim; }

static inline const struct var * design_tvar_item(const struct design *d, size_t k)
{
	assert(k <- design_tvar_count(d));
	return &d->tvars[k]->var;
}

const struct var *design_add_tvar(struct design *d, const char *name, const struct tvar_type *type, ...);

void design_get_tvar_matrix(const struct design *d, const double **x, const size_t **ind, size_t *nz);
const double *design_tvars(const struct design *d, size_t i);
const double *design_tvar(const struct design *d, const struct var *v, size_t i);




/* interactions */
const struct var *design_add_prod(struct design *d, const char *name, const struct var *u,
				  const struct var *v);

/* coefficients */
void coefs_init(struct coefs *c, const struct design *d);
void coefs_init_view(struct coefs *c, const struct design *d, const double *data);
void coefs_deinit(struct coefs *c);

/* algebra */
void design_mul(double alpha, const struct design *d,
		const struct coefs *c, double beta, double *y);
void design_traits_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y);
void design_tvars_mul(double alpha, const struct design *d,
		      const double *x, double beta, double *y);

void design_tmul(double alpha, const struct design *d, const double *x, double beta, struct coefs *c);
void design_traits_tmul(double alpha, const struct design *d, const double *x, double beta, double *y);
void design_tvars_tmul(double alpha, const struct design *d, const double *x, double beta, double *y);

void design_axpy(double alpha, const struct design *d, size_t i, struct coefs *c);
void design_traits_axpy(double alpha, const struct design *d, size_t i, double *y);
void design_tvars_axpy(double alpha, const struct design *d, size_t i, double *y);



/* internal functions (for use by tvar callbacks) */
double *design_make_active(struct design *d, struct tvar *v, size_t i);

static inline void design_update(struct design *d, struct tvar *v, size_t i, double t)
{
	assert(v->var.design == d);
	assert(i < design_count(d));
	assert(t >= d->tcur);

	if (t > deltaset_tlast(&d->deltaset, i)) {
		deltaset_update(&d->deltaset, i, t);

		if (t > deltaset_tlast(&v->deltaset, i))
			deltaset_update(&v->deltaset, i, t);
	}

}




#endif /* DESIGN_H */
