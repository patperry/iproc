#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include "memory.h"
#include "actors.h"

static int
compare_group_bucket (void *x1,
                      void *x2)
{
    iproc_group_bucket *bucket1 = x1;
    iproc_group_bucket *bucket2 = x2;
    size_t traits_hash1 = bucket1->traits_hash;
    size_t traits_hash2 = bucket2->traits_hash;

    if (traits_hash1 < traits_hash2) {
        return -1;
    } else if (traits_hash1 > traits_hash2) {
        return +1;
    } else {
        return 0;
    }
}

static int
compare_group (void *x1,
               void *x2)
{
    iproc_group *group1 = x1;
    iproc_group *group2 = x2;
    iproc_vector *traits1 = group1->traits;
    iproc_vector *traits2 = group2->traits;
    return iproc_vector_compare(traits1, traits2);
}

static int64_t
iproc_actors_insert_group (iproc_actors *actors,
                           iproc_vector *traits)
{
    assert(actors);
    assert(traits);

    /* first, find the right bucket */
    size_t traits_hash = iproc_vector_hash(traits);
    int64_t i = iproc_array_bsearch(actors->group_buckets, &traits_hash,
                                    compare_group_bucket);

    if (i < 0) { /* insert a new bucket if necessary */
        i = ~i;
        iproc_array *new_groups = iproc_array_new(sizeof(iproc_group));
        iproc_group_bucket new_bucket = { traits_hash, new_groups };
        iproc_array_insert(actors->group_buckets, i, &new_bucket);
    }
    
    iproc_group_bucket *bucket = &(iproc_array_index(actors->group_buckets,
                                                     iproc_group_bucket,
                                                     i));

    /* now, find the right group  */
    int64_t j = iproc_array_bsearch(bucket->groups, &traits,
                                    compare_group);

    if (j < 0) { /* insert a new group if necessary */
        j = ~j;
        iproc_vector *new_traits = iproc_vector_new_copy(traits);
        int64_t new_id = iproc_array_size(actors->group_traits);
        iproc_array_append(actors->group_traits, &new_traits);
        iproc_group new_group = { iproc_vector_ref(new_traits), new_id };

        iproc_array_insert(bucket->groups, j, &new_group);
    }
    
    /* lastly return the id */
    iproc_group *group = &(iproc_array_index(bucket->groups, iproc_group, j));
    return group->id;
}

static void
iproc_actors_group_ids_free (iproc_array *group_ids)
{
    iproc_array_unref(group_ids);
}

static void
iproc_actors_group_traits_free (iproc_array *group_traits)
{
    int64_t i, n;
    iproc_vector *x;

    if (group_traits) {
        n = iproc_array_size(group_traits);
        for (i = 0; i < n; i++) {
            x = iproc_array_index(group_traits, iproc_vector *, i);
            iproc_vector_unref(x);
        }
        iproc_array_unref(group_traits);
    }
}

static void
iproc_actors_group_bucket_deinit (iproc_group_bucket *bucket)
{
    if (!bucket)
        return;
    iproc_array *groups = bucket->groups;
    iproc_actors_group_traits_free(groups);
}

static void
iproc_actors_group_buckets_free (iproc_array *group_buckets)
{
    int64_t i, n;
    iproc_group_bucket *bucket;

    if (group_buckets) {
        n = iproc_array_size(group_buckets);
        for (i = 0; i < n; i++) {
            bucket = &(iproc_array_index(group_buckets, iproc_group_bucket, i));
            iproc_actors_group_bucket_deinit(bucket);
        }
        iproc_array_unref(group_buckets);
    }
}


static void
iproc_actors_free (iproc_actors *actors)
{
    if (actors) {
        iproc_actors_group_buckets_free(actors->group_buckets);
        iproc_actors_group_traits_free(actors->group_traits);
        iproc_actors_group_ids_free(actors->group_ids);
        iproc_free(actors);
    }
}



iproc_actors *
iproc_actors_new (int64_t      size,
                  iproc_vector *traits0)
{
    assert(0 <= size);
    assert(traits0);
    iproc_actors *actors = iproc_malloc(sizeof(*actors));
    if (!actors) return NULL;

    actors->group_ids = iproc_array_new(sizeof(int64_t));
    actors->group_traits = iproc_array_new(sizeof(iproc_vector *));
    actors->group_buckets = iproc_array_new(sizeof(iproc_group_bucket));
    iproc_refcount_init(&actors->refcount);

    if (!(actors->group_ids && actors->group_traits && actors->group_buckets)) {
        iproc_actors_free(actors);
        actors = NULL;
    }

    /* We rely on array_set_size clearing the tail to zeros. */
    iproc_array_set_size(actors->group_ids, size);

    iproc_actors_insert_group(actors, traits0);

    return actors;
}



iproc_actors *
iproc_actors_ref (iproc_actors *actors)
{
    if (actors) {
        iproc_refcount_get(&actors->refcount);
    }
    return actors;
}

static void
iproc_actors_release (iproc_refcount *refcount)
{
    iproc_actors *actors = container_of(refcount, iproc_actors, refcount);
    iproc_actors_free(actors);
}

void
iproc_actors_unref (iproc_actors *actors)
{
    if (!actors)
        return;

    iproc_refcount_put(&actors->refcount, iproc_actors_release);
}


int64_t
iproc_actors_size (iproc_actors *actors)
{
    assert(actors);
    return iproc_array_size(actors->group_ids);
}

int64_t
iproc_actors_dim (iproc_actors *actors)
{
    assert(actors);
    iproc_vector *x0 = iproc_actors_group_traits(actors, 0);
    int64_t n = iproc_vector_dim(x0);
    return n;
}

void
iproc_actors_set (iproc_actors *actors,
                  int64_t       actor_id,
                  iproc_vector *traits)
{
    assert(actors);
    assert(0 <= actor_id);
    assert(actor_id < iproc_actors_size(actors));
    assert(traits);
    assert(iproc_vector_dim(traits) == iproc_actors_dim(actors));

    int64_t group_id = iproc_actors_insert_group(actors, traits);
    iproc_array_index(actors->group_ids, int64_t, actor_id) = group_id;
}


iproc_vector *
iproc_actors_get (iproc_actors *actors,
                  int64_t       actor_id)
{
    assert(actors);
    assert(0 <= actor_id);
    assert(actor_id < iproc_actors_size(actors));

    int64_t g = iproc_actors_group(actors, actor_id);
    iproc_vector *x = iproc_actors_group_traits(actors, g);
    return x;
}

int64_t
iproc_actors_ngroup (iproc_actors *actors)
{
    assert(actors);
    return iproc_array_size(actors->group_traits);
}

int64_t
iproc_actors_group (iproc_actors *actors,
                    int64_t       actor_id)
{
    assert(actors);
    assert(0 <= actor_id);
    assert(actor_id < iproc_actors_size(actors));

    iproc_array *group_ids = actors->group_ids;
    int64_t g = iproc_array_index(group_ids, int64_t, actor_id);
    return g;
}


iproc_vector *
iproc_actors_group_traits (iproc_actors *actors,
                           int64_t       group_id)
{
    assert(actors);
    assert(0 <= group_id);
    assert(group_id < iproc_actors_ngroup(actors));

    iproc_vector *x = iproc_array_index(actors->group_traits,
                                        iproc_vector *,
                                        group_id);
    return x;
}


void
iproc_actors_mul (double        alpha,
                  iproc_trans   trans,
                  iproc_actors *actors,
                  iproc_vector *x,
                  double        beta,
                  iproc_vector *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_actors_dim(actors));

    int64_t n = iproc_actors_size(actors);
    int64_t i;
    iproc_vector *row;
    double dot, entry;

    if (beta == 0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1) {
        iproc_vector_scale(y, beta);
    }

    /* NOTE: this could be more efficient by using info about actor groups */
    if (trans == IPROC_TRANS_NOTRANS) {
        for (i = 0; i < n; i++) {
            row = iproc_actors_get(actors, i);
            dot = iproc_vector_dot(row, x);
            iproc_vector_inc(y, i, alpha * dot);
        }
    } else {
        for (i = 0; i < n; i++) {
            row = iproc_actors_get(actors, i);
            entry = iproc_vector_get(x, i);
            iproc_vector_acc(y, alpha * entry, row);
        }
    }
}


void
iproc_actors_muls (double         alpha,
                   iproc_trans    trans,
                   iproc_actors  *actors,
                   iproc_svector *x,
                   double         beta,
                   iproc_vector  *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_actors_dim(actors));

    int64_t n = iproc_actors_size(actors);
    int64_t i;
    iproc_vector *row;
    double dot, entry;

    if (beta == 0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1) {
        iproc_vector_scale(y, beta);
    }

    /* NOTE: this could be more efficient by using info about actor groups */
    if (trans == IPROC_TRANS_NOTRANS) {
        for (i = 0; i < n; i++) {
            row = iproc_actors_get(actors, i);
            dot = iproc_vector_sdot(row, x);
            iproc_vector_inc(y, i, alpha * dot);
        }
    } else {
        int64_t inz, nnz = iproc_svector_nnz(x);
        for (inz = 0; inz < nnz; inz++) {
            i = iproc_svector_nz(x, inz);
            row = iproc_actors_get(actors, i);
            entry = iproc_svector_nz_get(x, inz);
            iproc_vector_acc(y, alpha * entry, row);
        }
    }
}


void
iproc_actors_matmul (double        alpha,
                     iproc_trans   trans,
                     iproc_actors *actors,
                     iproc_matrix *x,
                     double        beta,
                     iproc_matrix *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(y) == iproc_actors_dim(actors));
    assert(iproc_matrix_ncol(x) == iproc_matrix_ncol(y));

    int64_t m = iproc_matrix_ncol(x);
    int64_t j;

    if (beta == 0) {
        iproc_matrix_set_all(y, 0.0);
    } else if (beta != 1) {
        iproc_matrix_scale(y, beta);
    }

    for (j = 0; j < m; j++) {
        iproc_vector_view xcol = iproc_matrix_col(x, j);
        iproc_vector_view ycol = iproc_matrix_col(y, j);
        iproc_actors_mul(alpha, trans, actors, &xcol.vector, 1.0, &ycol.vector);
    }
}