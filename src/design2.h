#ifndef DESIGN2_H
#define DESIGN2_H

#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include "deltaset.h"
#include "design.h"
#include "history.h"
#include "var.h"



struct var2 {
	struct var_meta meta;
	struct design2 *design;
	size_t index;
};

struct trait2 {
	struct var2 var;
	double *xi, *xj, *xij;
	size_t dimi, dimj;
};

struct tvar2 {
	struct var2 var;
	void *thunk;
	const struct tvar2_type *type;
	struct deltaset *deltaset;
	double *tcur;
};

struct tvar2_type {
	void (*init) (struct var_meta *meta, void **thunk, const char *name, struct design2 *d, va_list ap);
	void (*deinit) (struct var_meta *meta, void *thunk, struct design2 *d);
	double (*update) (struct tvar2 *tv, size_t i, double t0, const struct history *h);
};





struct design2 {
	struct history *history;
	size_t count1, count2;

	size_t *cohorts;
	size_t *cohort_reps;
	size_t ncohort;

	size_t trait_dim;
	double *trait_x;
	struct trait2 **traits;
	size_t ntrait;

	size_t tvar_dim;
	struct tvar2 **tvars;
	size_t ntvar, ntvar_max;

	size_t *tvar_ir;
	size_t *tvar_jc;
	double *tvar_x;
	size_t tvar_nnz, tvar_nnz_max;

	struct deltaset *deltaset;
	struct version_watch *history_version;
	double *tcur;
	double *tnext;
};


struct coefs2 {
	double *traits;
	double *tvars;
};


void design2_init(struct design2 *d, struct history *h, size_t count1, size_t count2);
void design2_deinit(struct design2 *d);



/* properties */
static inline struct history *design2_history(const struct design2 *d) { return d->history; }
static inline size_t design2_count1(const struct design2 *d) { return d->count1; }
static inline size_t design2_count2(const struct design2 *d) { return d->count2; }
static inline size_t design2_dim(const struct design2 *d) { return d->trait_dim + d->tvar_dim; }

const struct var2 *design2_var(const struct design2 *d, const char *name);



/* cohorts */
static inline size_t design2_cohort_count(const struct design2 *d) { return d->ncohort; }

static inline size_t design2_cohort(const struct design2 *d, size_t i)
{
	assert(i < design2_count1(d));
	return d->cohorts[i];
}

static inline size_t design2_cohort_rep(const struct design2 *d, size_t c)
{
	assert(c < design2_cohort_count(d));
	return d->cohort_reps[c];
}

static inline void design2_get_cohorts(const struct design2 *d, const size_t **cohorts,
				       const size_t **reps, size_t *ncohort)
{
	*cohorts = d->cohorts;
	*reps = d->cohort_reps;
	*ncohort = d->ncohort;
}



/* traits */
static inline size_t design2_trait_count(const struct design2 *d) { return d->ntrait; }
static inline size_t design2_trait_dim(const struct design2 *d) { return d->trait_dim; }

static inline const struct var2 *design2_trait_item(const struct design2 *d, size_t k)
{
	assert(k < design2_trait_count(d));
	return &d->traits[k]->var;
}


static inline const double *design2_trait_matrix(const struct design2 *d, size_t i)
{
	assert(i < design2_count1(d));
	size_t c = d->cohorts[i];
	size_t n = design2_count2(d);
	size_t dim = design2_trait_dim(d);
	const double *x = d->trait_x + c * n * dim;
	return x;
}

static inline const double *design2_traits(const struct design2 *d, size_t i, size_t j)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	const double *x = design2_trait_matrix(d, i);
	size_t dim = design2_trait_dim(d);
	return x + j * dim;
}

static inline const double *design2_trait(const struct design2 *d, const struct var2 *v, size_t i, size_t j)
{
	assert(v->design == d);
	assert(v->meta.type == VAR_TYPE_TRAIT);
	const double *x = design2_traits(d, i, j);
	return x + v->index;
}



/* tvars */
static inline size_t design2_tvar_count(const struct design2 *d) { return d->ntvar; }
static inline size_t design2_tvar_dim(const struct design2 *d) { return d->tvar_dim; }

static inline const struct var2 *design2_tvar_item(const struct design2 *d, size_t k)
{
	assert(k < design2_tvar_count(d));
	return &d->tvars[k]->var;
}

const struct var2 *design2_add_tvar(struct design2 *d, const char *name, const struct tvar2_type *type, ...);

void design2_get_tvar_matrix(const struct design2 *d, size_t i, const double **x, const size_t **ind, size_t *nz);
const double *design2_tvars(const struct design2 *d, size_t i, size_t j);
const double *design2_tvar(const struct design2 *d, const struct var2 *v, size_t i, size_t j);


double design2_next_time(const struct design2 *d, size_t i);
const struct deltaset *design2_changes(const struct design2 *d, size_t i);


/* interactions */
const struct var2 *design2_add_kron(struct design2 *d, const char *name, const struct var *i,
				    const struct var *j);
const struct var2 *design2_add_prod(struct design2 *d, const char *name, const struct var2 *u, const struct var2 *v);
//const struct var2 *design2_add_prod1(struct design2 *d, const char *name, const struct var *i, const struct var2 *v);
//const struct var2 *design2_add_prod2(struct design2 *d, const char *name, const struct var *j, const struct var2 *v);


/* algebra */
void design2_mul(double alpha, const struct design2 *d, size_t i, const struct coefs2 *c, double beta, double *y);
void design2_traits_mul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y);
void design2_tvars_mul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y);

void design2_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, struct coefs2 *c);
void design2_traits_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y);
void design2_tvars_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y);

void design2_axpy(double alpha, const struct design2 *d, size_t i, size_t j, struct coefs2 *c);
void design2_traits_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y);
void design2_tvars_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y);


/* internal functions (for use by tvar callbacks) */
double *design2_make_active(struct design2 *d, struct tvar2 *v, size_t i, size_t j);

static inline void design2_update(struct design2 *d, struct tvar2 *v, size_t i, size_t j, double t)
{
	assert(v->var.design == d);
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));
	assert(t <= d->tcur[i]);

	if (t > deltaset_tlast(&d->deltaset[i], j)) {
		deltaset_update(&d->deltaset[i], j, t);

		if (t > deltaset_tlast(&v->deltaset[i], j))
			deltaset_update(&v->deltaset[i], j, t);
	}

}



#endif /* DESIGN2_H */
