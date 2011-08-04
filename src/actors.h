#ifndef _ACTORS_H
#define _ACTORS_H

#include <stddef.h>

#include "array.h"
#include "intset.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"

struct dyad {
	ssize_t isend;
	ssize_t jrecv;
};

struct actor {
	ssize_t cohort;
};

struct cohort {
	struct array actors;
};

struct actors {
	struct array actors;
	struct array cohorts;

	/* deprecated */
	struct refcount refcount;
};

void actors_init(struct actors *a);
void actors_init_copy(struct actors *a, const struct actors *src);
void actors_deinit(struct actors *a);

void actors_clear(struct actors *a);

static inline ssize_t actors_count(const struct actors *a);
static inline const struct actor *actors_items(const struct actors *a);

static inline ssize_t actors_cohort_count(const struct actors *a);
static inline const struct cohort *actors_cohorts(const struct actors *a);

ssize_t actors_add(struct actors *a, ssize_t cohort);
ssize_t actors_add_cohort(struct actors *a);

void actors_mul(double alpha,
		enum trans_op trans,
		const struct actors *a,
		const struct vector *x, double beta, struct vector *y);

void actors_muls(double alpha,
		 enum trans_op trans,
		 const struct actors *a,
		 const struct svector *x, double beta, struct svector *y);

void actors_matmul(double alpha,
		   enum trans_op trans,
		   const struct actors *a,
		   const struct matrix *x, double beta, struct matrix *y);

/* deprecated */
struct actors *actors_alloc(void);
struct actors *actors_ref(struct actors *a);
void actors_free(struct actors *a);

/* inline function definitions */
ssize_t actors_count(const struct actors *a)
{
	assert(a);
	return array_count(&a->actors);
}

ssize_t actors_cohort_count(const struct actors *a)
{
	assert(a);
	return array_count(&a->cohorts);
}

const struct actor *actors_items(const struct actors *a)
{
	assert(a);
	return array_to_ptr(&a->actors);
}

const struct cohort *actors_cohorts(const struct actors *a)
{
	assert(a);
	return array_to_ptr(&a->cohorts);
}

#endif /* _ACTORS_H */
