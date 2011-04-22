#ifndef _IPROC_ACTORS_H
#define _IPROC_ACTORS_H

#include <stddef.h>

#include "darray.h"
#include "refcount.h"

#include "hashset.h"
#include "matrix.h"
#include "svector.h"
#include "vector.h"

#include "cohort.h"

typedef struct actors iproc_actors;

/*
 
  
 

struct _iproc_actors {
    struct hashset  cohorts;
    struct darray   actors;
    ssize_t         dim; 
    struct refcount refcount;
};
*/

typedef struct _iproc_group iproc_group;
struct _iproc_group {
	struct vector *traits;	/* traits must be the first member */
	int64_t id;
};

typedef struct _iproc_group_bucket iproc_group_bucket;
struct _iproc_group_bucket {
	size_t traits_hash;	/* traits_hash must be the first member */
	struct darray groups;
};

struct actor {
	struct cohort *cohort;
};

struct actors {
	ssize_t dim;
	struct darray actors;
	struct hashset cohorts;

	/* deprecated */
	struct darray group_ids;
	struct darray group_traits;
	struct darray group_buckets;
	struct refcount refcount;
};

/* makes a copy of traits0 */
struct actors *actors_alloc(ssize_t dim);
struct actors *actors_ref(struct actors *a);
void actors_free(struct actors *a);

ssize_t actors_size(const struct actors *a);
ssize_t actors_cohorts_size(const struct actors *a);
ssize_t actors_dim(const struct actors *a);

const struct actor *actors_at(const struct actors *a, ssize_t actor_id);

/* makes a copy of traits */
ssize_t actors_add(iproc_actors * actors, const struct vector *traits);
const struct vector *actors_traits(const struct actors *a, ssize_t actor_id);

void actors_mul(double alpha,
		iproc_trans trans,
		const iproc_actors * actors,
		const struct vector *x, double beta, struct vector *y);
void actors_muls(double alpha,
		 iproc_trans trans,
		 const iproc_actors * actors,
		 const iproc_svector * x, double beta, struct vector *y);

void actors_matmul(double alpha,
		   iproc_trans trans,
		   const iproc_actors * actors,
		   const iproc_matrix * x, double beta, iproc_matrix * y);

/* deprecated */
int64_t iproc_actors_group(const iproc_actors * actors, int64_t actor_id);
const struct vector *iproc_actors_group_traits(const iproc_actors * actors,
					       int64_t group_id);

#endif /* _IPROC_ACTORS_H */
