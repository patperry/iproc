#ifndef CONSTR_H
#define CONSTR_H

#include "blas.h"
#include <assert.h>
#include <stddef.h>


struct constr {
	size_t dim;
	double *wts;
	double *vals;
	size_t n, nmax;
};

/* init/deinit */
void constr_init(struct constr *c, size_t dim);
void constr_init_copy(struct constr *c, const struct constr *c0);
void constr_deinit(struct constr *c);

/* properties */
static inline size_t constr_dim(const struct constr *c);
static inline size_t constr_count(const struct constr *c);
static inline const double *constr_wts(const struct constr *c, size_t k);
static inline double constr_val(const struct constr *c, size_t k);
static inline const double *constr_all_wts(const struct constr *c);
static inline const double *constr_all_vals(const struct constr *c);

/* adding constraints */
int constr_add(struct constr *c, const double *wts, double val);
int constr_add_set(struct constr *c, size_t i, double val);
int constr_add_eq(struct constr *c, size_t i1, size_t i2);
size_t constr_add_identify(struct constr *c, const double *hess,
			   enum blas_uplo uplo);


/* inline function definitions */
size_t constr_dim(const struct constr *c)
{
	return c->dim;
}

size_t constr_count(const struct constr *c)
{
	return c->n;
}

const double *constr_wts(const struct constr *c, size_t k)
{
	assert(k < constr_count(c));
	const double *wts = constr_all_wts(c);
	size_t dim = constr_dim(c);
	return wts + k * dim;
}

double constr_val(const struct constr *c, size_t k)
{
	assert(k < constr_count(c));
	const double *vals = constr_all_vals(c);
	return vals[k];
}

const double *constr_all_wts(const struct constr *c)
{
	return c->wts;
}

const double *constr_all_vals(const struct constr *c)
{
	return c->vals;
}



#endif /* CONSTR_H */