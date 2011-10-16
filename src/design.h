#ifndef DESIGN_H
#define DESIGN_H

#include "blas.h"
#include "messages.h"
#include "svector.h"
#include "vector.h"


struct var_type;		// forward declaration
struct design_var;

struct design {
	size_t count;
	double *intvls;
	size_t nintvl;

	size_t dim;
	int has_effects;

	size_t trait_off;
	size_t trait_dim;
	double *traits;
	char **trait_names;

	size_t dvar_off;
	size_t dvar_dim;
	struct design_var *dvars;
	size_t ndvar, ndvar_max;
};

void design_init(struct design *d, size_t count,
		 const double *intvls, size_t nintvl,
		 const double *traits, size_t trait_dim,
		 const char *const *trait_names);
void design_deinit(struct design *d);

static inline size_t design_count(const struct design *d);
static inline const double *design_intervals(const struct design *d);
static inline size_t design_interval_count(const struct design *d);

static inline size_t design_dim(const struct design *d);

static inline int design_has_effects(const struct design *d);
void design_set_has_effects(struct design *d, int has_effects);
static inline size_t design_effects_index(const struct design *d);

static inline size_t design_traits_index(const struct design *d);
static inline size_t design_traits_dim(const struct design *d);
static inline const double *design_traits(const struct design *d);
static inline const char *const *design_trait_names(const struct design *d);


void design_add_dvar(struct design *d, const struct var_type *type,
		     void *params);
ptrdiff_t design_dvar_index(const struct design *d,
			    const struct var_type *type);
static inline void design_get_dvars(const struct design *d,
				    const struct design_var **dvarsp,
				    size_t *np);

static inline size_t design_dvars_index(const struct design *d);
static inline size_t design_dvars_dim(const struct design *d);

void design_mul0(double alpha,
		      enum blas_trans trans,
		      const struct design *d,
		      const struct vector *x, double beta, struct vector *y);
void design_muls0(double alpha,
		       enum blas_trans trans,
		       const struct design *d,
		       const struct svector *x, double beta, struct vector *y);

/* inline funciton definitions */
size_t design_count(const struct design *d)
{
	assert(d);
	return d->count;
}

const double *design_intervals(const struct design *d)
{
	assert(d);
	return d->intvls;
}

size_t design_interval_count(const struct design *d)
{
	assert(d);
	return d->nintvl;
}

size_t design_dim(const struct design *d)
{
	return d->dim;
}

int design_has_effects(const struct design *d)
{
	return d->has_effects;
}

size_t design_effects_index(const struct design *d)
{
	assert(design_has_effects(d));
	(void)d;
	return 0;
}

size_t design_traits_index(const struct design *d)
{
	return d->trait_off;
}

size_t design_traits_dim(const struct design *d)
{
	return d->trait_dim;
}

const double *design_traits(const struct design *d)
{
	return d->traits;
}

const char * const *design_trait_names(const struct design *d)
{
	return (const char * const *)d->trait_names;
}

void design_get_dvars(const struct design *d,
				    const struct design_var **dvarsp,
				    size_t *np)
{
	*dvarsp = d->dvars;
	*np = d->ndvar;
}

size_t design_dvars_index(const struct design *d)
{
	return d->dvar_off;
}

size_t design_dvars_dim(const struct design *d)
{
	return d->dvar_dim;
}

#endif /* DESIGN_H */
