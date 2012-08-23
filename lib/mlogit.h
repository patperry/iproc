#ifndef MLOGIT_H
#define MLOGIT_H

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include "blas.h"
#include "catdist.h"

#define MLOGIT_COV_UPLO	BLAS_UPPER

struct mlogit {
	struct catdist values;
	double *x;
	double *beta;
	double *offset;
	double *mean;
	double *cov;		// true covariance is cov / exp(log_cov_scale)
	double log_cov_scale;

	double *mean_diff;
	double *cov_diff;
	double *cat_buf;
	double *dim_buf1;
	double *dim_buf2;

	size_t dim;
	double mean_err;
	double cov_err;
	double log_cov_scale_err;
};

void mlogit_init(struct mlogit *m, size_t ncat, size_t dim);
void mlogit_deinit(struct mlogit *m);

static inline size_t mlogit_ncat(const struct mlogit *m);
static inline size_t mlogit_dim(const struct mlogit *m);

static inline double *mlogit_coefs(const struct mlogit *m);
static inline double *mlogit_offset(const struct mlogit *m);
static inline double *mlogit_x(const struct mlogit *m);
static inline double *mlogit_mean(const struct mlogit *m);
static inline double *mlogit_cov(const struct mlogit *m, double *cov_scale);
static inline struct catdist *mlogit_values(const struct mlogit *m);

void mlogit_set_coefs(struct mlogit *m, const double *beta);
void mlogit_set_all_offset(struct mlogit *m, const double *offset);
void mlogit_set_offset(struct mlogit *m, size_t i, double offset);
void mlogit_set_all_x(struct mlogit *m, const double *x);
void mlogit_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		  const double *dx, size_t ndx);

int _mlogit_check(const struct mlogit *m);

size_t mlogit_ncat(const struct mlogit *m)
{
	return catdist_ncat(&m->values);
}

size_t mlogit_dim(const struct mlogit *m)
{
	return m->dim;
}

double *mlogit_coefs(const struct mlogit *m)
{
	return m->beta;
}

double *mlogit_offset(const struct mlogit *m)
{
	return m->offset;
}

double *mlogit_x(const struct mlogit *m)
{
	return m->x;
}

double *mlogit_mean(const struct mlogit *m)
{
	return m->mean;
}

double *mlogit_cov(const struct mlogit *m, double *cov_scale)
{
	*cov_scale = exp(m->log_cov_scale);
	return m->cov;
}

struct catdist *mlogit_values(const struct mlogit *m)
{
	return &((struct mlogit *)m)->values;
}

#endif /* MLOGIT_H */
