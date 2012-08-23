#ifndef MLOGITAUG_H
#define MLOGITAUG_H

#include <assert.h>
#include <stddef.h>
#include "catdist.h"
#include "catdist1.h"
#include "mlogit.h"

struct mlogitaug {
	const struct mlogit *base;
	struct catdist1 dist;

	double *beta;
	size_t dim;

	size_t *ind;
	double *x1;
	double *offset;
	size_t nz, nzmax;
};


void mlogitaug_init(struct mlogitaug *m1, const struct mlogitaug *base, size_t dim);
void mlogitaug_deinit(struct mlogitaug *m);
void mlogitaug_clear(struct mlogitaug *m);

static struct mlogitaug *mlogitaug_base(const struct mlogitaug *m1);
static inline size_t mlogitaug_ncat(const struct mlogitaug *m1);
static inline size_t mlogitaug_dim(const struct mlogitaug *m1);

double mlogitaug_offset(const struct mlogitaug *m1);
double mlogitaug_doffset(const struct mlogitaug *m1, size_t i);
void mlogitaug_set_doffset(struct mlogitaug *m1, size_t i, double doffset);

double *mlogitaug_x1(const struct mlogitaug *m1, size_t i);
double *mlogitaug_mean0(const struct mlogitaug *m);
double *mlogitaug_mean1(const struct mlogitaug *m);
double *mlogitaug_cov0(const struct mlogitaug *m, double *cov_scale);
double *mlogitaug_cov1(const struct mlogitaug *m, double *cov_scale);

struct mlogitaug *mlogitaug_values(const struct mlogitaug *m1);

void mlogitaug_inc_dx(struct mlogitaug *m1, size_t i,
			   const size_t *jdx, const double *dx, size_t ndx);

int _mlogitaug_check(const struct mlogitaug *m1);

#endif /* MLOGITAUG_H */
