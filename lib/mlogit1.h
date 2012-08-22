#ifndef MLOGIT1_H
#define MLOGIT1_H

#include <assert.h>
#include <stdlib.h>
#include "mlogit.h"

struct mlogit1 {
	const struct mlogit *parent;

	double *deta;
	size_t *ind;
	size_t nz, nzmax;
};

void mlogit1_init(struct mlogit1 *m1, const struct mlogit *parent);
void mlogit1_deinit(struct mlogit1 *m1);
void mlogit1_clear(struct mlogit1 *m1);

static inline size_t mlogit1_ncat(const struct mlogit1 *m1);

double mlogit1_prob(const struct mlogit1 *m1, size_t i);
double mlogit1_lprob(const struct mlogit1 *m1, size_t i);

double mlogit1_eta(const struct mlogit1 *m1, size_t i);
double mlogit_deta(const struct mlogit1 *m1, size_t i);
void mlogit1_set_deta(struct mlogit1 *m1, size_t i, double deta);

double mlogit1_psi(const struct mlogit1 *m1);
double mlogit_dpsi(const struct mlogit1 *m1);
static inline void mlogit1_get_deta(const struct mlogit1 *m1,
				    const double **deta, const size_t **ind,
				    size_t *nz);

int _mlogit1_check(const struct mlogit1 *m);

size_t mlogit1_ncat(const struct mlogit1 *m1)
{
	return mlogit_ncat(m1->parent);
}

void mlogit1_get_deta(const struct mlogit1 *m1, const double **deta,
		      const size_t **ind, size_t *nz)
{
	*deta = m1->deta;
	*ind = m1->ind;
	*nz = m1->nz;
}

#endif /* MLOGIT1_H */
