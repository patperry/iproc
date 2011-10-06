#ifndef STRATA_H
#define STRATA_H

#include <stddef.h>
#include "hashset.h"

struct strata {
	size_t dim;
	size_t nlevel;
	size_t nlevel_max;
	double **levels;
	struct hashset levels_set;
};

void strata_init(struct strata *s, size_t dim);
void strata_deinit(struct strata *s);
void strata_clear(struct strata *s);

static inline size_t strata_dim(const struct strata *s);
static inline size_t strata_count(const struct strata *s);

size_t strata_add(struct strata *s, const double *level);
ptrdiff_t strata_find(const struct strata *s, const double *level);
static inline double **strata_levels(const struct strata *s);

/* inline function definitions */
size_t strata_dim(const struct strata *s)
{
	return s->dim;
}

size_t strata_count(const struct strata *s)
{
	return s->nlevel;
}

double **strata_levels(const struct strata *s)
{
	return s->levels;
}

#endif /* STRATA_H */
