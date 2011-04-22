#include "port.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "compare.h"
#include "actors.h"

DEFINE_COMPARE_FN(int64_compare, int64_t)
DEFINE_COMPARE_FN(size_compare, size_t)

#define group_bucket_compare size_compare

static bool iproc_actors_reserve(iproc_actors * actors, ssize_t size);

static int group_compare(const void *x, const void *y)
{
	return vector_compare(((iproc_group *) x)->traits,
			      ((iproc_group *) y)->traits);
}

static ssize_t
iproc_actors_insert_group(iproc_actors * actors, const struct vector *traits)
{
	assert(actors);
	assert(traits);

	/* first, find the right bucket */
	size_t traits_hash = vector_hash(traits);
	ssize_t i = darray_binary_search(&actors->group_buckets,
					 &traits_hash,
					 group_bucket_compare);

	if (i < 0) {		/* insert a new bucket if necessary */
		i = ~i;
		iproc_group_bucket new_bucket;
		new_bucket.traits_hash = traits_hash;
		darray_init(&new_bucket.groups, sizeof(iproc_group));
		darray_insert(&actors->group_buckets, i, &new_bucket);
	}

	iproc_group_bucket *bucket = darray_at(&actors->group_buckets, i);

	/* now, find the right group  */
	ssize_t j = darray_binary_search(&bucket->groups,
					 &traits,
					 group_compare);

	if (j < 0) {		/* insert a new group if necessary */
		j = ~j;
		struct vector *new_traits = vector_new_copy(traits);
		int64_t new_id = darray_size(&actors->group_traits);
		darray_push_back(&actors->group_traits, &new_traits);
		iproc_group new_group = { new_traits, new_id };

		darray_insert(&bucket->groups, j, &new_group);
	}

	/* lastly return the id */
	iproc_group *group = darray_at(&bucket->groups, j);
	return group->id;
}

static void iproc_actors_group_deinit(iproc_group * group)
{
	vector_free(group->traits);
}

static void iproc_actors_group_bucket_deinit(iproc_group_bucket * bucket)
{
	ssize_t i, n;
	struct darray *groups = &bucket->groups;
	iproc_group *group;

	n = darray_size(groups);
	for (i = 0; i < n; i++) {
		group = darray_at(groups, i);
		iproc_actors_group_deinit(group);
	}

	darray_deinit(groups);
}

static void iproc_actors_group_buckets_deinit(struct darray *group_buckets)
{
	ssize_t i, n;
	iproc_group_bucket *bucket;

	n = darray_size(group_buckets);
	for (i = 0; i < n; i++) {
		bucket = darray_at(group_buckets, i);
		iproc_actors_group_bucket_deinit(bucket);
	}
	darray_deinit(group_buckets);
}

static uint32_t cohortp_traits_hash(const void *cohortp)
{
	return vector_hash(cohort_traits(*(struct cohort **)cohortp));
}

static bool cohortp_traits_equals(const void *cohortp1, const void *cohortp2)
{
	return vector_equals(cohort_traits(*(struct cohort **)cohortp1),
			     cohort_traits(*(struct cohort **)cohortp2));
}

static bool actors_init(struct actors *actors, ssize_t dim)
{
	assert(actors);
	assert(dim >= 0);

	bool ok;

	ok = darray_init(&actors->actors, sizeof(struct actor));
	if (!ok)
		goto fail_actors;

	ok = hashset_init(&actors->cohorts, cohortp_traits_hash,
			  cohortp_traits_equals, sizeof(struct cohort *));
	if (!ok)
		goto fail_cohorts;

	/* old */
	ok = darray_init(&actors->group_ids, sizeof(int64_t));
	if (!ok)
		goto fail_group_ids;

	ok = darray_init(&actors->group_traits, sizeof(struct vector *));
	if (!ok)
		goto fail_group_traits;

	ok = darray_init(&actors->group_buckets, sizeof(iproc_group_bucket));
	if (!ok)
		goto fail_group_buckets;

	ok = refcount_init(&actors->refcount);
	if (!ok)
		goto fail_refcount;
	/* end of old */
	
	actors->dim = dim;
	return true;

	refcount_deinit(&actors->refcount);
fail_refcount:
	darray_deinit(&actors->group_buckets);
fail_group_buckets:
	darray_deinit(&actors->group_traits);
fail_group_traits:
	darray_deinit(&actors->group_ids);
fail_group_ids:
	hashset_deinit(&actors->cohorts);
fail_cohorts:
	darray_deinit(&actors->actors);	
fail_actors:
	return false;
}

static void cohorts_clear(struct hashset *cohorts)
{
	struct hashset_iter it;
	hashset_iter_init(cohorts, &it);
	while (hashset_iter_advance(cohorts, &it)) {
		struct cohort **cp = hashset_iter_current(cohorts, &it);
		cohort_free(*cp);
	}
	hashset_iter_deinit(cohorts, &it);
}

static void cohorts_deinit(struct hashset *cohorts)
{
	cohorts_clear(cohorts);
	hashset_deinit(cohorts);
}

static void actors_deinit(struct actors *actors)
{
	assert(actors);

	cohorts_deinit(&actors->cohorts);
	darray_deinit(&actors->actors);

	/* old */
	refcount_deinit(&actors->refcount);
	iproc_actors_group_buckets_deinit(&actors->group_buckets);
	darray_deinit(&actors->group_traits);
	darray_deinit(&actors->group_ids);
}

struct actors *actors_alloc(ssize_t dim)
{
	assert(dim >= 0);
	struct actors *actors;

	if ((actors = malloc(sizeof(*actors)))) {
		if (actors_init(actors, dim))
			return actors;
		free(actors);
	}
	return NULL;
}

void actors_free(iproc_actors * actors)
{
	if (!actors)
		return;

	if (refcount_put(&actors->refcount, NULL)) {
		actors_deinit(actors);
		free(actors);
	}
}

struct actors *actors_ref(struct actors *actors)
{
	assert(actors);
	refcount_get(&actors->refcount);
	return actors;
}

ssize_t actors_size(const struct actors *a)
{
	assert(a);
	return darray_size(&a->actors);
}

ssize_t actors_cohorts_size(const struct actors *a)
{
	assert(a);
	return hashset_size(&a->cohorts);
}
 
ssize_t actors_dim(const struct actors *actors)
{
	assert(actors);
	return actors->dim;
}

const struct actor *actors_at(const struct actors *a, ssize_t actor_id)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_size(a));
	return darray_at(&a->actors, actor_id);
}

const struct cohort *actors_cohort(const struct actors *a, ssize_t actor_id)
{
	return actors_at(a, actor_id)->cohort;
}

/* deprecated */
bool iproc_actors_reserve(iproc_actors * actors, ssize_t size)
{
	assert(actors);
	assert(size >= 0);

	return darray_reserve(&actors->group_ids, size);
}

/* deprecated */
ssize_t actors_add_old(iproc_actors * actors, const struct vector *traits)
{
	assert(actors);
	assert(vector_size(traits) == actors_dim(actors));

	ssize_t id = -1;

	if (iproc_actors_reserve(actors, actors_size(actors) + 1)) {
		id = actors_size(actors);

		int64_t group_id = iproc_actors_insert_group(actors, traits);

		darray_push_back(&actors->group_ids, &group_id);
	}

	return id;
}


static struct cohort * actors_get_cohort(struct actors *actors,
					 const struct vector *traits)
{
	struct cohort *key = container_of(traits, struct cohort, traits);
	struct hashset_pos pos;
	struct cohort *new_cohort;
	struct cohort **cohortp;
	
	if ((cohortp = hashset_find(&actors->cohorts, &key, &pos)))
		return *cohortp;
	
	if ((new_cohort = cohort_new(traits))) {
		if (hashset_insert(&actors->cohorts, &pos, &new_cohort)) {
			return new_cohort;
		}
		cohort_free(new_cohort);
	}
	return NULL;
}

ssize_t actors_add(struct actors *actors, const struct vector *traits)
{
	assert(actors);
	assert(vector_size(traits) == actors_dim(actors));

	struct actor a;
	ssize_t id;
	bool ok;

	a.cohort = actors_get_cohort(actors, traits);
	if (!a.cohort)
		goto fail_get_cohort;

	ok = darray_push_back(&actors->actors, &a);
	if (!ok)
		goto fail_push_back;
	
	id = darray_size(&actors->actors) - 1;
	ok = cohort_add(a.cohort, id);
	if (!ok)
		goto fail_cohort_add;
	
	goto out;
	
fail_cohort_add:
	darray_pop_back(&actors->actors);
fail_push_back:
fail_get_cohort:
	id = -1;
out:
	actors_add_old(actors, traits);
	return id;
}

const struct vector *actors_traits(const iproc_actors * a,
				   ssize_t actor_id)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_size(a));

	const struct actor *actor = actors_at(a, actor_id);
	return cohort_traits(actor->cohort);
}

/* deprecated
int64_t actors_cohort(const iproc_actors * actors, int64_t actor_id)
{
	assert(actors);
	assert(0 <= actor_id);
	assert(actor_id < actors_size(actors));

	const struct darray *group_ids = &actors->group_ids;
	int64_t g = *(int64_t *)darray_at(group_ids, actor_id);
	return g;
}
 */

/* deprecated */
const struct vector *iproc_actors_group_traits(const iproc_actors * actors,
					       int64_t group_id)
{
	assert(actors);
	assert(0 <= group_id);
	assert(group_id < actors_cohorts_size(actors));

	struct vector *x =
	    *(struct vector **)darray_at(&actors->group_traits, group_id);
	return x;
}

void
actors_mul(double alpha,
	   iproc_trans trans,
	   const struct actors *a,
	   const struct vector *x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_size(x) == actors_dim(a));
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_size(y) == actors_size(a));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_size(x) == actors_size(a));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_size(y) == actors_dim(a));

	const struct vector *row;
	double alpha_dot, scale;
	struct hashset_iter it;
	struct cohort *c;
	struct cohort_iter c_it;
	ssize_t id;
	
	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	hashset_iter_init(&a->cohorts, &it);
	if (trans == IPROC_TRANS_NOTRANS) {
		while (hashset_iter_advance(&a->cohorts, &it)) {
			c = *(struct cohort **)hashset_iter_current(&a->cohorts, &it);
			row = cohort_traits(c);
			alpha_dot = alpha * vector_dot(row, x);
			
			cohort_iter_init(c, &c_it);
			while (cohort_iter_advance(c, &c_it)) {
				id = cohort_iter_current(c, &c_it);
				*vector_at(y, id) += alpha_dot;
			}
			cohort_iter_deinit(c, &c_it);
		}
	} else {
		while (hashset_iter_advance(&a->cohorts, &it)) {
			c = *(struct cohort **)hashset_iter_current(&a->cohorts, &it);
			row = cohort_traits(c);
			scale = 0.0;
			
			cohort_iter_init(c, &c_it);
			while (cohort_iter_advance(c, &c_it)) {
				id = cohort_iter_current(c, &c_it);
				scale += *vector_at(x, id);
			}
			cohort_iter_deinit(c, &c_it);
			
			vector_axpy(alpha * scale, row, y);
		}

		
	}
	hashset_iter_deinit(&a->cohorts, &it);
}

void
actors_muls(double alpha,
	    iproc_trans trans,
	    const struct actors *a,
	    const iproc_svector * x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || iproc_svector_dim(x) == actors_dim(a));
	assert(trans != IPROC_TRANS_NOTRANS
	       || vector_size(y) == actors_size(a));
	assert(trans == IPROC_TRANS_NOTRANS
	       || iproc_svector_dim(x) == actors_size(a));
	assert(trans == IPROC_TRANS_NOTRANS
	       || vector_size(y) == actors_dim(a));

	const struct vector *row;
	double alpha_dot, entry;
	struct hashset_iter it;
	struct cohort *c;
	struct cohort_iter c_it;
	ssize_t id;

	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	if (trans == IPROC_TRANS_NOTRANS) {
		hashset_iter_init(&a->cohorts, &it);		
		while (hashset_iter_advance(&a->cohorts, &it)) {
			c = *(struct cohort **)hashset_iter_current(&a->cohorts, &it);
			row = cohort_traits(c);
			alpha_dot = alpha * vector_dots(row, x);
			
			cohort_iter_init(c, &c_it);
			while (cohort_iter_advance(c, &c_it)) {
				id = cohort_iter_current(c, &c_it);
				*vector_at(y, id) += alpha_dot;
			}
			cohort_iter_deinit(c, &c_it);
		}
		hashset_iter_deinit(&a->cohorts, &it);		
	} else {
		/* NOTE: this could potentially be made more effecient by
		 * using the cohort structure.  Below, we assume that
		 * the sparsity in x is more important.
		 */
		ssize_t inz, nnz = iproc_svector_nnz(x);
		for (inz = 0; inz < nnz; inz++) {
			id = iproc_svector_nz(x, inz);
			row = actors_traits(a, id);
			entry = iproc_svector_nz_get(x, inz);
			vector_axpy(alpha * entry, row, y);
		}
	}	
}

void
actors_matmul(double alpha,
	      iproc_trans trans,
	      const struct actors *a,
	      const iproc_matrix * x, double beta, iproc_matrix * y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != IPROC_TRANS_NOTRANS
	       || iproc_matrix_nrow(x) == actors_dim(a));
	assert(trans != IPROC_TRANS_NOTRANS
	       || iproc_matrix_nrow(y) == actors_size(a));
	assert(trans == IPROC_TRANS_NOTRANS
	       || iproc_matrix_nrow(x) == actors_size(a));
	assert(trans == IPROC_TRANS_NOTRANS
	       || iproc_matrix_nrow(y) == actors_dim(a));
	assert(iproc_matrix_ncol(x) == iproc_matrix_ncol(y));

	int64_t m = iproc_matrix_ncol(x);
	int64_t j;

	if (beta == 0) {
		iproc_matrix_set_all(y, 0.0);
	} else if (beta != 1) {
		iproc_matrix_scale(y, beta);
	}

	struct vector xcol;
	struct vector ycol;

	for (j = 0; j < m; j++) {
		vector_init_matrix_col(&xcol, x, j);
		vector_init_matrix_col(&ycol, y, j);
		actors_mul(alpha, trans, a, &xcol, 1.0, &ycol);
	}
}
