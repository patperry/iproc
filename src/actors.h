#ifndef _ACTORS_H
#define _ACTORS_H

#include <stddef.h>

#include "darray.h"
#include "hashset.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"

#include "cohort.h"

struct actor {
	struct cohort *cohort;
};

struct actors {
	ssize_t dim;
	struct darray actors;
	struct hashset cohorts;
	struct refcount refcount;
};

struct actors *actors_alloc(ssize_t dim);
struct actors *actors_ref(struct actors *a);
void actors_free(struct actors *a);
bool actors_init(struct actors *actors, ssize_t dim);
void actors_deinit(struct actors *a);

void actors_clear(struct actors *a);

ssize_t actors_size(const struct actors *a);
ssize_t actors_cohorts_size(const struct actors *a);
ssize_t actors_dim(const struct actors *a);

/* makes a copy of traits */
ssize_t actors_add(struct actors *a, const struct vector *traits);

const struct actor *actors_at(const struct actors *a, ssize_t actor_id);
const struct vector *actors_traits(const struct actors *a, ssize_t actor_id);
const struct cohort *actors_cohort(const struct actors *a, ssize_t actor_id);

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

#endif /* _ACTORS_H */
