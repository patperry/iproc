#ifndef _STRATA_H
#define _STRATA_H

#include "array.h"
#include "intset.h"
#include "intmap.h"
#include "vector.h"


struct stratum {
	struct vector traits;
	struct intset ids;
};

struct strata {
	ssize_t dim;
	struct array items;
	struct intmap trait_hashes;
};

void strata_init(struct strata *s, ssize_t dim);
void strata_deinit(struct strata *s);

static inline ssize_t strata_dim(const struct strata *s);
static inline ssize_t strata_count(const struct strata *s);
static inline struct stratum *strata_items(const struct strata *s);

ssize_t strata_add(struct strata *s, intptr_t id, const struct vector *traits);

static inline struct vector *stratum_traits(const struct stratum *s);
static inline ssize_t stratum_count(const struct stratum *s);
static inline intptr_t *stratum_items(const struct stratum *s);


/* inline function definitions */
ssize_t strata_dim(const struct strata *s)
{
	assert(s);
	return s->dim;
}

ssize_t strata_count(const struct strata *s)
{
	assert(s);
	return array_count(&s->items);
}

struct stratum *strata_items(const struct strata *s)
{
	assert(s);
	return array_to_ptr(&s->items);
}

struct vector *stratum_traits(const struct stratum *s)
{
	assert(s);
	return &((struct stratum *)s)->traits;
}

ssize_t stratum_count(const struct stratum *s)
{
	assert(s);
	return intset_count(&s->ids);
}

intptr_t *stratum_items(const struct stratum *s)
{
	assert(s);
	return intset_items(&s->ids);
}



#endif /* _STRATA_H */
