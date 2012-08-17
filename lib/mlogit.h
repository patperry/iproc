#ifndef MLOGIT_H
#define MLOGIT_H

#include <assert.h>
#include <stddef.h>


struct mlogit {
	size_t ncat;

	double *eta; /* linear predictors */
	size_t *eta_order;
	size_t *eta_rank;

	double eta_max;
	double eta_tail; // \sum_{ i != i_max } exp{ eta(i) - eta_max }
	double eta_tail_err;
	double psi_shift; /* shifted CGF: log[ \sum_i exp{ eta(i) - eta_max } ] */
	double tol;
};


void mlogit_init(struct mlogit *m, size_t ncat);
void mlogit_deinit(struct mlogit *m);
void mlogit_clear(struct mlogit *m);

static inline size_t mlogit_ncat(const struct mlogit *m);

static inline double mlogit_eta(const struct mlogit *m, size_t i);
static inline double mlogit_prob(const struct mlogit *m, size_t i);
static inline double mlogit_lprob(const struct mlogit *m, size_t i);
static inline double mlogit_psi(const struct mlogit *m);

void mlogit_set_eta(struct mlogit *m, size_t i, double eta);
void mlogit_set_all_eta(struct mlogit *m, const double *eta);


void _mlogit_check_invariants(const struct mlogit *m);

size_t mlogit_ncat(const struct mlogit *m)
{
	return m->ncat;
}

double mlogit_eta(const struct mlogit *m, size_t i)
{
	assert(i < mlogit_ncat(m));
	return m->eta[i];
}

double mlogit_prob(const struct mlogit *m, size_t i)
{
	assert(i < mlogit_ncat(m));
	return exp(mlogit_lprob(m, i));
}

double mlogit_lprob(const struct mlogit *m, size_t i)
{
	assert(i < mlogit_ncat(m));
	return (m->eta[i] - m->eta_max) - m->psi_shift;
}

double mlogit_psi(const struct mlogit *m)
{
	assert(!isnan(m->eta_max + m->psi_shift));
	return m->eta_max + m->psi_shift;
}



#endif /* MLOGIT_H */
