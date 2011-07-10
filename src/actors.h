#ifndef _ACTORS_H
#define _ACTORS_H

#include <stddef.h>

#include "array.h"
#include "hashset.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"

#include "cohort.h"

struct actor {
	ssize_t cohort;
};

struct actors {
	ssize_t dim;
	struct array actors;
	struct array cohorts;
	struct intmap trait_hashes; /* (hash(x_i), { i_1, ... i_k }) */
	
	/* deprecated */
	struct refcount refcount;
};

struct dyad {
	ssize_t isend;
	ssize_t jrecv;
};

struct actors *actors_alloc(ssize_t dim);
struct actors *actors_ref(struct actors *a);
void actors_free(struct actors *a);
void actors_init(struct actors *actors, ssize_t dim);
void actors_init_copy(struct actors *actors, const struct actors *src);
void actors_init_matrix(struct actors *actors, enum trans_op trans, const struct matrix *matrix);	// copies x
void actors_deinit(struct actors *a);

void actors_clear(struct actors *a);

static inline ssize_t actors_count(const struct actors *a);
static inline ssize_t actors_cohort_count(const struct actors *a);
static inline ssize_t actors_dim(const struct actors *a);

/* makes a copy of traits */
void actors_add(struct actors *a, const struct vector *traits);

static inline const struct actor *actors_item(const struct actors *a, ssize_t actor_id);
static inline const struct vector *actors_traits(const struct actors *a, ssize_t actor_id);
static inline const struct cohort *actors_cohort(const struct actors *a, ssize_t cohort_id);

void actors_mul(double alpha,
		enum trans_op trans,
		const struct actors *a,
		const struct vector *x, double beta, struct vector *y);
void actors_muls(double alpha,
		 enum trans_op trans,
		 const struct actors *a,
		 const struct svector *x, double beta, struct vector *y);

void actors_matmul(double alpha,
		   enum trans_op trans,
		   const struct actors *a,
		   const struct matrix *x, double beta, struct matrix *y);

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

ssize_t actors_dim(const struct actors *actors)
{
	assert(actors);
	return actors->dim;
}

const struct actor *actors_item(const struct actors *a, ssize_t actor_id)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_count(a));
	return array_item(&a->actors, actor_id);
	//const struct actor *base = array_to_ptr(&a->actors);
	//return base + actor_id;
}

const struct vector *actors_traits(const struct actors *a, ssize_t actor_id)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_count(a));
	
	const struct actor *actor = actors_item(a, actor_id);
	const struct cohort *cohort = actors_cohort(a, actor->cohort);
	return cohort_traits(cohort);
}

const struct cohort *actors_cohort(const struct actors *a, ssize_t cohort_id)
{
	assert(a);
	assert(0 <= cohort_id && cohort_id < actors_cohort_count(a));
	return array_item(&a->cohorts, cohort_id);
	//const struct cohort *base = array_to_ptr(&a->cohorts);
	//return base + cohort_id;
}



#endif /* _ACTORS_H */
