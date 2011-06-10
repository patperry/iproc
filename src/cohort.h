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

struct cohort_iter {
	struct intset_iter members_it;
};

#define COHORT_KEY(it) ((ssize_t)INTSET_KEY(it.members_it))
#define COHORT_FOREACH(it, c) \
	for ((it) = cohort_iter_make(c); cohort_iter_advance(&(it));)


/* create, destroy */
struct cohort *cohort_alloc(const struct vector *x);	// copies x
void cohort_free(struct cohort *c);
void cohort_init(struct cohort *c, const struct vector *x);	// copies x
void cohort_deinit(struct cohort *c);

/* traits */
static inline ssize_t cohort_dim(const struct cohort *c);
const struct vector *cohort_traits(const struct cohort *c);

/* size/query */
static inline ssize_t cohort_count(const struct cohort *c);
bool cohort_contains(const struct cohort *c, ssize_t id);

/* modification */
bool cohort_add(struct cohort *c, ssize_t id);
bool cohort_remove(struct cohort *c, ssize_t id);

/* iteration */
struct cohort_iter cohort_iter_make(const struct cohort *c);
bool cohort_iter_advance(struct cohort_iter *it);
void cohort_iter_reset(struct cohort_iter *it);


/* inline funcion definitions */
ssize_t cohort_dim(const struct cohort *c)
{
	assert(c);
	return vector_dim(cohort_traits(c));
}

ssize_t cohort_count(const struct cohort *c)
{
	assert(c);
	return intset_count(&c->members);
}


#endif /* _COHORT */
