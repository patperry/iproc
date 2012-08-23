#include "port.h"
#include "mlogitaug.h"

void mlogitaug_init(struct mlogitaug *m1, const struct mlogitaug *base,
		size_t dim)
{
}

void mlogitaug_deinit(struct mlogitaug *m)
{
}

double *mlogitaug_coefs(const struct mlogitaug *m1)
{
	return m1->beta;
}

double mlogitaug_offset(const struct mlogitaug *m1, size_t i)
{
	return 0.0;
}

double *mlogitaug_x(const struct mlogitaug *m1, size_t i)
{
	return NULL;
}

void mlogitaug_set_coefs(struct mlogitaug *m1, const double *beta)
{
}

void mlogitaug_set_offset(struct mlogitaug *m1, size_t i, double offset)
{
}

void mlogitaug_inc_x(struct mlogitaug *m, size_t i, const size_t *jdx,
		const double *dx, size_t ndx)
{
}

struct catdist1 *mlogitaug_dist(const struct mlogitaug *m1)
{
	return NULL;
}

int mlogitaug_check(const struct mlogitaug *m1)
{
	return 0;
}
