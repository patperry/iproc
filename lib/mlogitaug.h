#ifndef MLOGITAUG_H
#define MLOGITAUG_H

#include <assert.h>
#include <stddef.h>
#include "catdist.h"
#include "catdist1.h"
#include "mlogit.h"

struct mlogitaug_work {
	double *cov_full;
	double *xbuf;
};

struct mlogitaug {
	const struct mlogit *base;
	struct catdist1 dist;

	double *beta;
	size_t dim;

	size_t *ind;
	double *offset;
	double *x;
	double *deta;
	size_t nz, nzmax;

	double *mean;
	double *cov;

	double *base_mean;
	double *base_dmean;
	double *base_cov;

	double *cross_cov;

	struct mlogitaug_work *work;
	int free_work;
};

void mlogitaug_work_init(struct mlogitaug_work *work, size_t base_dim,
			 size_t aug_dim);
void mlogitaug_work_deinit(struct mlogitaug_work *work);


void mlogitaug_init(struct mlogitaug *m1, const struct mlogit *base,
		    size_t dim, struct mlogitaug_work *work);
void mlogitaug_deinit(struct mlogitaug *m1);

static inline size_t mlogitaug_ncat(const struct mlogitaug *m1);
static inline size_t mlogitaug_dim(const struct mlogitaug *m1);
static inline const struct mlogit *mlogitaug_base(const struct mlogitaug *m1);

double *mlogitaug_coefs(const struct mlogitaug *m1);
double mlogitaug_offset(const struct mlogitaug *m1, size_t i);
double *mlogitaug_x(const struct mlogitaug *m1, size_t i);

void mlogitaug_set_coefs(struct mlogitaug *m1, const double *beta);
void mlogitaug_set_offset(struct mlogitaug *m1, size_t i, double offset);
void mlogitaug_set_all_offset(struct mlogitaug *m1, const size_t *i, const double *offset, size_t nz);
void mlogitaug_inc_x(struct mlogitaug *m1, size_t i, const size_t *jdx,
		const double *dx, size_t ndx);
void mlogitaug_set_x(struct mlogitaug *m1, size_t i, const double *x);
void mlogitaug_set_all_x(struct mlogitaug *m1, const size_t *i, const double *x, size_t nz);

struct catdist1 *mlogitaug_dist(const struct mlogitaug *m1);

double *mlogitaug_mean(const struct mlogitaug *m1);
double *mlogitaug_base_mean(const struct mlogitaug *m1);


double *mlogitaug_cov(const struct mlogitaug *m1);
double *mlogitaug_base_cov(const struct mlogitaug *m1);
double *mlogitaug_cross_cov(const struct mlogitaug *m1);


void mlogitaug_update_cache(struct mlogitaug *m1);
struct catdist1 *mlogitaug_cached_dist(const struct mlogitaug *m1);
double *mlogitaug_cached_mean(const struct mlogitaug *m1);
double *mlogitaug_cached_base_mean(const struct mlogitaug *m1);
double *mlogitaug_cached_cov(const struct mlogitaug *m1);
double *mlogitaug_cached_base_cov(const struct mlogitaug *m1);
double *mlogitaug_cached_cross_cov(const struct mlogitaug *m1);

int mlogitaug_check(const struct mlogitaug *m1);



size_t mlogitaug_ncat(const struct mlogitaug *m1)
{
	return mlogit_ncat(m1->base);
}

size_t mlogitaug_dim(const struct mlogitaug *m1)
{
	return m1->dim;
}

const struct mlogit *mlogitaug_base(const struct mlogitaug *m1)
{
	return m1->base;
}


#endif /* MLOGITAUG_H */
