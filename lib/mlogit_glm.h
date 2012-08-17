#ifndef MLOGIT_GLM_H
#define MLOGIT_GLM_H

#include <assert.h>
#include <stddef.h>

#include "mlogit.h"


struct mlogit_glm {
	struct mlogit mlogit;
	struct dmatrix x;
	double *beta;
	double *mean;
	double *cov;
	size_t dim;
};


void mlogit_glm_init(struct mlogit_glm *m, size_t ncat, size_t dim);
void mlogit_glm_deinit(struct mlogit_glm *m);

void mlogit_glm_set_coefs(struct mlogit_glm *m, const double *beta);
void mlogit_glm_set_x(struct mlogit_glm *m, const double *x);
void mlogit_glm_inc_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *idx, size_t ndx);


struct mlogit *mlogit_glm_values(const struct mlogit_glm *m);
double *mlogit_glm_mean(const struct mlogit_glm *m);
double *mlogit_glm_cov(const struct mlogit_glm *m);



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
 
 void mlogit_update(struct mlogit *m, size_t i, double deta);
 
 void mlogit_mean_init(struct mlogit_mean *m, size_t dim, const double *mean0);
 void mlogit_mean_deinit(struct mlogit_mean *m);
 
 void mlogit_mean_update(struct mlogit_mean *m, const struct mlogit *mlogit,
 const double *x1, const double *dx, const struct vpattern *ix);
 
 void mlogit_cov_update(struct mlogit_cov *m, const struct mlogit *mlogit,
 const struct mlogit_mean *mean);
 
 */


#endif /* MLOGIT_GLM_H */