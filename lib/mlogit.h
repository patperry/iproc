#ifndef MLOGIT_H
#define MLOGIT_H

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include "blas.h"
#include "catdist.h"

#define MLOGIT_COV_UPLO	BLAS_UPPER

struct mlogit_work {
	double *mean_diff;
	double *cov_diff;
	double *cov_diff_full;
	double *cat_buf;
	double *dim_buf1;
	double *dim_buf2;
};

struct mlogit {
	struct catdist dist;
	double *x;
	double *beta;
	double *offset;
	double *mean;
	double *cov;		// true covariance is cov / exp(log_cov_scale)
	double log_cov_scale;

	size_t dim;
	double mean_err;
	double cov_err;
	double log_cov_scale_err;

	int moments;

	struct mlogit_work *work;
	int free_work;

};

void mlogit_work_init(struct mlogit_work *work, size_t ncat, size_t dim);
void mlogit_work_deinit(struct mlogit_work *work);


void mlogit_init(struct mlogit *m, size_t ncat, size_t dim, struct mlogit_work *work);
void mlogit_deinit(struct mlogit *m);

static inline size_t mlogit_ncat(const struct mlogit *m);
static inline size_t mlogit_dim(const struct mlogit *m);

static inline int mlogit_moments(const struct mlogit *m);
void mlogit_set_moments(struct mlogit *m, int k);

static inline double *mlogit_coefs(const struct mlogit *m);
static inline double *mlogit_offset(const struct mlogit *m);
static inline double *mlogit_x(const struct mlogit *m);
static inline double *mlogit_mean(const struct mlogit *m);
static inline double *mlogit_cov(const struct mlogit *m, double *cov_scale);
static inline struct catdist *mlogit_dist(const struct mlogit *m);

void mlogit_set_coefs(struct mlogit *m, const double *beta);
void mlogit_set_all_offset(struct mlogit *m, const double *offset);
void mlogit_set_offset(struct mlogit *m, size_t i, double offset);
void mlogit_set_all_x(struct mlogit *m, const double *x);
void mlogit_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		  const double *dx, size_t ndx);

int mlogit_check(const struct mlogit *m);

size_t mlogit_ncat(const struct mlogit *m)
{
	return catdist_ncat(&m->dist);
}

size_t mlogit_dim(const struct mlogit *m)
{
	return m->dim;
}

int mlogit_moments(const struct mlogit *m)
{
	return m->moments;
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
	assert(m->moments >= 1);
	return m->mean;
}

double *mlogit_cov(const struct mlogit *m, double *cov_scale)
{
	assert(m->moments >= 2);
	*cov_scale = exp(m->log_cov_scale);
	return m->cov;
}

struct catdist *mlogit_dist(const struct mlogit *m)
{
	return &((struct mlogit *)m)->dist;
}

#endif /* MLOGIT_H */
