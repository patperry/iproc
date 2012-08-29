#ifndef CATDIST1_H
#define CATDIST1_H

#include <assert.h>
#include <stdlib.h>
#include "catdist.h"

struct catdist1 {
	const struct catdist *parent;

	size_t *ind;
	double *deta;
	size_t nz, nzmax;

	double cached_dpsi;
};

void catdist1_init(struct catdist1 *c1, const struct catdist *parent);
void catdist1_deinit(struct catdist1 *c1);

static inline size_t catdist1_ncat(const struct catdist1 *c1);
static inline const struct catdist *catdist1_parent(const struct catdist1 *c1);

double catdist1_prob(const struct catdist1 *c1, size_t i);
double catdist1_lprob(const struct catdist1 *c1, size_t i);

double catdist1_eta(const struct catdist1 *c1, size_t i);
double catdist1_deta(const struct catdist1 *c1, size_t i);
void catdist1_set_deta(struct catdist1 *c1, size_t i, double deta);
void catdist1_set_all_deta(struct catdist1 *c1, const size_t *ind, const double *deta, size_t nz);
static inline void catdist1_get_deta(const struct catdist1 *c1,
				    const size_t **ind, const double **deta,
				    size_t *nz);

double catdist1_psi(const struct catdist1 *c1);
double catdist1_dpsi(const struct catdist1 *c1);

int catdist1_check(const struct catdist1 *c1);


double catdist1_psi(const struct catdist1 *c1);
double catdist1_dpsi(const struct catdist1 *c1);

/* fast (unsafe) operations; uses cached value of dpsi */
void catdist1_update_cache(struct catdist1 *c1);
double catdist1_cached_prob(const struct catdist1 *c1, size_t i);
double catdist1_cached_lprob(const struct catdist1 *c1, size_t i);
double catdist1_cached_psi(const struct catdist1 *c1);
double catdist1_cached_dpsi(const struct catdist1 *c1);




size_t catdist1_ncat(const struct catdist1 *c1)
{
	return catdist_ncat(c1->parent);
}

const struct catdist *catdist1_parent(const struct catdist1 *c1)
{
	return c1->parent;
}

void catdist1_get_deta(const struct catdist1 *c1,
		      const size_t **ind,
		      const double **deta,
		      size_t *nz)
{
	*ind = c1->ind;
	*deta = c1->deta;
	*nz = c1->nz;
}

#endif /* CATDIST1_H */
