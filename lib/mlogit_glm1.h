#ifndef MLOGIT_GLM_AUG_H
#define MLOGIT_GLM_AUG_H

#include <assert.h>
#include <stddef.h>
#include "mlogit.h"
#include "mlogit1.h"
#include "mlogit_glm.h"

struct mlogit_glm_aug {
	const struct mlogit_glm *base;
	struct mlogit1 values;

	double *dx;
	double *doffset;
	size_t *ind;
	size_t nz, nzmax;
};

void mlogit_glm_aug_init(struct mlogit_glm_aug *m1,
			 const struct mlogit_glm *base);
void mlogit_glm_aug_deinit(struct mlogit_glm_aug *m);
void mlogit_glm_aug_clear(struct mlogit_glm_aug *m);

static inline size_t mlogit_glm_aug_ncat(const struct mlogit_glm_aug *m1);
static inline size_t mlogit_glm_aug_dim(const struct mlogit_glm_aug *m1);

double mlogit_glm_aug_offset(const struct mlogit_glm_aug *m1, size_t i);
double mlogit_glm_aug_doffset(const struct mlogit_glm_aug *m1, size_t i);

double *mlogit_glm_aug_dx(const struct mlogit_glm_aug *m1, size_t i);
double *mlogit_glm_mean(const struct mlogit_glm *m);
double *mlogit_glm_cov(const struct mlogit_glm *m, double *cov_scale);
struct mlogit1 *mlogit_glm_aug_values(const struct mlogit_glm_aug *m1);

void mlogit_glm_aug_set_doffset(struct mlogit_glm_aug *m1, size_t i,
				double doffset);
void mlogit_glm_aug_inc_dx(struct mlogit_glm_aug *m1, size_t i,
			   const size_t *jdx, const double *dx, size_t ndx);

int _mlogit_glm_aug_check(const struct mlogit_glm_aug *m1);

#endif /* MLOGIT_GLM_AUG_H */
