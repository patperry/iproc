#ifndef MLOGIT_H
#define MLOGIT_H

#include "blas.h"
#include "sblas.h"

struct mlogit {
	size_t n;
	
	double *eta; /* linear predictors */
	size_t *eta_order;
	size_t *eta_rank;

	double eta_max;
	double eta_tail; // \sum_{ i != i_max } exp{ eta(i) - eta_max }
	double phi; /* shifted CGF: log[ \sum_i exp{ eta(i) - eta_max } ] */

	// last update
	double eta0;
	double deta;
	double expm1_deta;
};

struct mlogit_mean {
	size_t dim;
	double *mean; /* expected covariates */
	double *xbuf;
};

struct mlogit_cov {
	double *cov; /* covariance of covariates (packed) */
	enum blas_uplo uplo;
};

void mlogit_init(struct mlogit *m, size_t n, const double *eta0);
void mlogit_deinit(struct mlogit *m);

void mlogit_update(struct mlogit *m, size_t i, double deta);

void mlogit_mean_init(struct mlogit_mean *m, size_t dim, const double *mean0);
void mlogit_mean_deinit(struct mlogit_mean *m);

void mlogit_mean_update(struct mlogit_mean *m, const struct mlogit *mlogit,
			const double *x1, const double *dx, const struct vpattern *ix);

void mlogit_cov_update(struct mlogit_cov *m, const struct mlogit *mlogit,
		       const struct mlogit_mean *mean);


#endif /* MLOGIT_H */
