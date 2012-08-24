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

void catdist_init(struct catdist *c, size_t ncat);
void catdist_deinit(struct catdist *c);

static inline size_t catdist_ncat(const struct catdist *c);

static inline double catdist_eta(const struct catdist *c, size_t i);
static inline double catdist_prob(const struct catdist *c, size_t i);
static inline double catdist_lprob(const struct catdist *c, size_t i);
static inline double catdist_psi(const struct catdist *c);

void catdist_set_eta(struct catdist *c, size_t i, double eta);
void catdist_set_all_eta(struct catdist *c, const double *eta);

int catdist_check(const struct catdist *c);

size_t catdist_ncat(const struct catdist *c)
{
	return c->ncat;
}

double catdist_eta(const struct catdist *c, size_t i)
{
	assert(i < catdist_ncat(c));
	return c->eta[i];
}

double catdist_prob(const struct catdist *c, size_t i)
{
	assert(i < catdist_ncat(c));
	return exp(catdist_lprob(c, i));
}

double catdist_lprob(const struct catdist *c, size_t i)
{
	assert(i < catdist_ncat(c));
	return (c->eta[i] - c->eta_max) - c->psi_shift;
}

double catdist_psi(const struct catdist *c)
{
	return c->eta_max + c->psi_shift;
}

#endif /* CATDIST_H */
