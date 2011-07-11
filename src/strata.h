#ifndef _STRATA_H
#define _STRATA_H

#include "array.h"
#include "intset.h"
#include "intmap.h"
#include "vector.h"


struct strata {
	ssize_t dim;
	struct array levels;
	struct intmap level_hashes;
};

void strata_init(struct strata *s, ssize_t dim);
void strata_deinit(struct strata *s);

static inline ssize_t strata_dim(const struct strata *s);
static inline ssize_t strata_count(const struct strata *s);

ssize_t strata_add(struct strata *s, const struct vector *level);
ssize_t strata_lookup(const struct strata *s, const struct vector *level);
static inline struct vector *strata_levels(const struct strata *s);

/* inline function definitions */
ssize_t strata_dim(const struct strata *s)
{
	assert(s);
	return s->dim;
}

ssize_t strata_count(const struct strata *s)
{
	assert(s);
	return array_count(&s->levels);
}

struct vector *strata_levels(const struct strata *s)
{
	assert(s);
	return array_to_ptr(&s->levels);
}


#endif /* _STRATA_H */
