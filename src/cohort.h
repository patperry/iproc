#ifndef _COHORT
#define _COHORT

#include "darray.h"
#include "vector.h"

/* A set of actors with the same traits
 */

struct cohort {
	struct vector traits;
	struct darray members;
};

struct cohort_iter {
	ssize_t index;
};

/* create, destroy */
struct cohort *cohort_new(const struct vector *x);	// copies x
void cohort_free(struct cohort *c);
bool cohort_init(struct cohort *c, const struct vector *x);	// copies x
void cohort_deinit(struct cohort *c);

/* traits */
ssize_t cohort_dim(const struct cohort *c);
const struct vector *cohort_traits(const struct cohort *c);

/* size/query */
bool cohort_empty(const struct cohort *c);
ssize_t cohort_size(const struct cohort *c);
bool cohort_contains(const struct cohort *c, ssize_t id);

/* modification */
bool cohort_add(struct cohort *c, ssize_t id);
void cohort_remove(struct cohort *c, ssize_t id);

/* iteration */
void cohort_iter_init(const struct cohort *c, struct cohort_iter *it);
void cohort_iter_deinit(const struct cohort *c, struct cohort_iter *it);
bool cohort_iter_advance(const struct cohort *c, struct cohort_iter *it);
ssize_t cohort_iter_current(const struct cohort *c,
			    const struct cohort_iter *it);

#endif /* _COHORT */
