#ifndef CATDIST_H
#define CATDIST_H

#include <assert.h>		// assert
#include <math.h>		// exp, isnan
#include <stddef.h>		// size_t

/* A categorical distribution, using the natural exponential family
 * parametrization.
 */
struct catdist {
	size_t ncat;

	double *eta;		/* linear predictors */
	size_t *eta_order;
	size_t *eta_rank;

	double eta_max;
	double eta_tail;	// \sum_{ i != i_max } exp{ eta(i) - eta_max }
	double eta_tail_err;
	double psi_shift;	/* shifted CGF: log[ \sum_i exp{ eta(i) - eta_max } ] */
	double tol;
};

void catdist_init(struct catdist *m, size_t ncat);
void catdist_deinit(struct catdist *m);

static inline size_t catdist_ncat(const struct catdist *m);

static inline double catdist_eta(const struct catdist *m, size_t i);
static inline double catdist_prob(const struct catdist *m, size_t i);
static inline double catdist_lprob(const struct catdist *m, size_t i);
static inline double catdist_psi(const struct catdist *m);

void catdist_set_eta(struct catdist *m, size_t i, double eta);
void catdist_set_all_eta(struct catdist *m, const double *eta);

int _catdist_check(const struct catdist *m);

size_t catdist_ncat(const struct catdist *m)
{
	return m->ncat;
}

double catdist_eta(const struct catdist *m, size_t i)
{
	assert(i < catdist_ncat(m));
	return m->eta[i];
}

double catdist_prob(const struct catdist *m, size_t i)
{
	assert(i < catdist_ncat(m));
	return exp(catdist_lprob(m, i));
}

double catdist_lprob(const struct catdist *m, size_t i)
{
	assert(i < catdist_ncat(m));
	return (m->eta[i] - m->eta_max) - m->psi_shift;
}

double catdist_psi(const struct catdist *m)
{
	return m->eta_max + m->psi_shift;
}

#endif /* CATDIST_H */
