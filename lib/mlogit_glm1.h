#ifndef MLOGIT_GLM1_H
#define MLOGIT_GLM1_H

#include <assert.h>
#include <stddef.h>
#include "mlogit.h"
#include "mlogit1.h"
#include "mlogit_glm.h"




struct mlogit_glm1 {
	const struct mlogit_glm *parent;
	struct mlogit1 values;

	double *dx;
	double *doffset;
	size_t *ind;
	size_t nz, nzmax;
};


void mlogit_glm1_init(struct mlogit_glm1 *m1, const struct mlogit_glm *parent);
void mlogit_glm1_deinit(struct mlogit_glm1 *m1);
void mlogit_glm1_clear(struct mlogit_glm1 *m1);

static inline size_t mlogit_glm1_ncat(const struct mlogit_glm1 *m1);
static inline size_t mlogit_glm1_dim(const struct mlogit_glm1 *m1);

double mlogit_glm1_offset(const struct mlogit_glm1 *m1, size_t i);
double mlogit_glm1_doffset(const struct mlogit_glm1 *m1, size_t i);

double *mlogit_glm1_dx(const struct mlogit_glm1 *m1, size_t i);
double *mlogit_glm_mean(const struct mlogit_glm *m);
double *mlogit_glm_cov(const struct mlogit_glm *m, double *cov_scale);
struct mlogit1 *mlogit_glm1_values(const struct mlogit_glm1 *m1);

void mlogit_glm1_set_doffset(struct mlogit_glm1 *m1, size_t i, double doffset);
void mlogit_glm1_inc_dx(struct mlogit_glm1 *m1, size_t i, const size_t *jdx, const double *dx, size_t ndx);


int _mlogit_glm1_check(const struct mlogit_glm1 *m1);


#endif /* MLOGIT_GLM1_H */
