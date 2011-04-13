#ifndef _COHORT
#define _COHORT

#include "intset.h"
#include "vector.h"

/* A set of actors with the same traits
 */

struct cohort {
	struct vector traits;
	struct intset members;
};

/* create, destroy */
struct cohort *cohort_new(const struct vector *x);	// copies x
void cohort_free(struct cohort *c);
bool cohort_init(struct cohort *c, const struct vector *x);	// copies x
void cohort_deinit(struct cohort *c);

/* traits */
static inline ssize_t cohort_dim(const struct cohort *c);
static inline const struct vector *cohort_traits(const struct cohort *c);

/* size/query */
static inline bool cohort_empty(const struct cohort *c);
static inline ssize_t cohort_size(const struct cohort *c);
static inline ssize_t cohort_at(const struct cohort *c, ssize_t i);
static inline bool cohort_contains(const struct cohort *c, ssize_t id);

/* modification */
bool cohort_add(struct cohort *c, ssize_t id);
void cohort_remove(struct cohort *c, ssize_t id);

/* inline function definitions */
ssize_t cohort_dim(const struct cohort *c)
{
	assert(c);
	return vector_size(cohort_traits(c));
}

const struct vector *cohort_traits(const struct cohort *c)
{
	assert(c);
	return &c->traits;
}

bool cohort_empty(const struct cohort * c)
{
	assert(c);
	return intset_empty(&c->members);
}

ssize_t cohort_size(const struct cohort * c)
{
	assert(c);
	return intset_size(&c->members);
}

bool cohort_contains(const struct cohort * c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	return intset_contains(&c->members, id);
}

ssize_t cohort_at(const struct cohort * c, ssize_t i)
{
	assert(c);
	assert(0 <= i && i < cohort_size(c));
	return intset_at(&c->members, i);
}

#endif /* _COHORT */
