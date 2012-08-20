#ifndef MLOGIT_GLM_H
#define MLOGIT_GLM_H

#include <assert.h>
#include <stddef.h>
#include "blas.h"
#include "mlogit.h"


#define MLOGIT_GLM_COV_UPLO	BLAS_LOWER


struct mlogit_glm {
	struct mlogit values;
	double *x;
	double *beta;
	double *mean;
	double *cov;
	double *cat_buf;
	double *dim_buf;
	size_t dim;
	double mean_err;
};


void mlogit_glm_init(struct mlogit_glm *m, size_t ncat, size_t dim);
void mlogit_glm_deinit(struct mlogit_glm *m);

static inline size_t mlogit_glm_ncat(const struct mlogit_glm *m);
static inline size_t mlogit_glm_dim(const struct mlogit_glm *m);
static inline double *mlogit_glm_coefs(const struct mlogit_glm *m);
static inline double *mlogit_glm_x(const struct mlogit_glm *m);
static inline double *mlogit_glm_mean(const struct mlogit_glm *m);
static inline double *mlogit_glm_cov(const struct mlogit_glm *m);
static inline struct mlogit *mlogit_glm_values(const struct mlogit_glm *m);

void mlogit_glm_set_coefs(struct mlogit_glm *m, const double *beta);
void mlogit_glm_set_all_x(struct mlogit_glm *m, const double *x);
void mlogit_glm_inc_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx);



void _mlogit_glm_check_invariants(const struct mlogit_glm *m);


size_t mlogit_glm_ncat(const struct mlogit_glm *m)
{
	return mlogit_ncat(&m->values);
}

size_t mlogit_glm_dim(const struct mlogit_glm *m)
{
	return m->dim;
}

double *mlogit_glm_coefs(const struct mlogit_glm *m)
{
	return m->beta;
}

double *mlogit_glm_x(const struct mlogit_glm *m)
{
	return m->x;
}

double *mlogit_glm_mean(const struct mlogit_glm *m)
{
	return m->mean;
}

double *mlogit_glm_cov(const struct mlogit_glm *m)
{
	return m->cov;
}

struct mlogit *mlogit_glm_values(const struct mlogit_glm *m)
{
	return &((struct mlogit_glm *)m)->values;
}




/*
 
 struct mlogit_mean {
 size_t dim;
 double *mean; // expected covariates
 double *xbuf;
 };
 
 struct mlogit_cov {
 double *cov; // covariance of covariates (packed)
 enum blas_uplo uplo;
 };
 
 
 */


#endif /* MLOGIT_GLM_H */