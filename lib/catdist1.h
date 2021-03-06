#ifndef CATDIST1_H
#define CATDIST1_H

#include <assert.h>
#include <stdlib.h>
#include "uintset.h"
#include "catdist.h"


struct catdist1 {
	struct catdist *parent;
	struct version_watch parent_watch;

	struct uintset ind;
	double *deta;
	size_t nzmax;

	double dpsi;
	int cached;
};

void catdist1_init(struct catdist1 *c1, struct catdist *parent);
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
	uintset_get_vals(&c1->ind, ind, nz);
	*deta = c1->deta;
}

#endif /* CATDIST1_H */
