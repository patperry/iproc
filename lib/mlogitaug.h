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
	double *off;
	double *x;
	size_t nz, nzmax;
};


void mlogitaug_init(struct mlogitaug *m1, const struct mlogitaug *base,
		size_t dim);
void mlogitaug_deinit(struct mlogitaug *m);

static inline size_t mlogitaug_ncat(const struct mlogitaug *m1);
static inline size_t mlogitaug_dim(const struct mlogitaug *m1);
static inline struct mlogit *mlogitaug_base(const struct mlogitaug *m1);

double *mlogitaug_coefs(const struct mlogitaug *m1);
double mlogitaug_offset(const struct mlogitaug *m1);
double *mlogitaug_x(const struct mlogitaug *m1, size_t i);

void mlogitaug_set_coefs(struct mlogitaug *m1, const double *beta);
void mlogitaug_set_offset(struct mlogitaug *m1, size_t i, double offset);
void mlogitaug_inc_x(struct mlogitaug *m, size_t i, const size_t *jdx,
		const double *dx, size_t ndx);

// double *mlogitaug_mean(const struct mlogitaug *m);
// double *mlogitaug_base_mean(const struct mlogitaug *m);
//
// double *mlogitaug_cov(const struct mlogitaug *m, double *scale);
// double *mlogitaug_base_cov(const struct mlogitaug *m, double *scale);
// double *mlogitaug_cross_cov(const struct mlogitaug *m, double *scale);

struct catdist1 *mlogitaug_dist(const struct mlogitaug *m1);

int mlogitaug_check(const struct mlogitaug *m1);

#endif /* MLOGITAUG_H */
