#include "port.h"

#include <assert.h>
#include "compare.h"
#include "memory.h"
#include "actors.h"

#define group_bucket_compare size_compare

static bool iproc_actors_reserve (iproc_actors *actors, ssize_t size);


static int group_compare (const void *x, const void *y)
{
    return vector_compare(((iproc_group *)x)->traits,
                          ((iproc_group *)y)->traits);
}


static ssize_t
iproc_actors_insert_group (iproc_actors *actors, const struct vector *traits)
{
    assert(actors);
    assert(traits);

    /* first, find the right bucket */
    size_t traits_hash = vector_hash(traits);
    ssize_t i = darray_binary_search(&actors->group_buckets,
                                     &traits_hash,
                                     group_bucket_compare);

    if (i < 0) { /* insert a new bucket if necessary */
        i = ~i;
        iproc_group_bucket new_bucket;
        new_bucket.traits_hash = traits_hash;
        darray_init(&new_bucket.groups, iproc_group);
        darray_insert(&actors->group_buckets, i, &new_bucket);
    }
    
    iproc_group_bucket *bucket = &darray_index(&actors->group_buckets,
                                               iproc_group_bucket,
                                               i);

    /* now, find the right group  */
    ssize_t j = darray_binary_search(&bucket->groups,
                                     &traits,
                                     group_compare);

    if (j < 0) { /* insert a new group if necessary */
        j = ~j;
        struct vector *new_traits = vector_new_copy(traits);
        int64_t new_id = darray_size(&actors->group_traits);
        darray_push_back(&actors->group_traits, &new_traits);
        iproc_group new_group = { new_traits, new_id };

        darray_insert(&bucket->groups, j, &new_group);
    }
    
    /* lastly return the id */
    iproc_group *group = &darray_index(&bucket->groups, iproc_group, j);
    return group->id;
}


static void
iproc_actors_group_deinit (iproc_group *group)
{
    vector_free(group->traits);
}


static void
iproc_actors_group_bucket_deinit (iproc_group_bucket *bucket)
{
    ssize_t i, n;
    struct darray *groups = &bucket->groups;
    iproc_group *group;
    
    n = darray_size(groups);
    for (i = 0; i < n; i++) {
        group = &darray_index(groups, iproc_group, i);
        iproc_actors_group_deinit(group);
    }
    
    darray_deinit(groups);
}


static void
iproc_actors_group_buckets_deinit (struct darray *group_buckets)
{
    ssize_t i, n;
    iproc_group_bucket *bucket;

    n = darray_size(group_buckets);
    for (i = 0; i < n; i++) {
        bucket = &darray_index(group_buckets, iproc_group_bucket, i);
        iproc_actors_group_bucket_deinit(bucket);
    }
    darray_deinit(group_buckets);
}


static void
iproc_actors_deinit (iproc_actors *actors)
{
    iproc_actors_group_buckets_deinit(&actors->group_buckets);
    darray_deinit(&actors->group_traits);
    darray_deinit(&actors->group_ids);
}


static void
iproc_actors_free (iproc_actors *actors)
{
    if (actors) {
        iproc_actors_deinit(actors);
        iproc_free(actors);
    }
}



iproc_actors *
iproc_actors_new (ssize_t dim)
{
    assert(dim >= 0);
    iproc_actors *actors = malloc(sizeof(*actors));

    if (actors) {
        if (darray_init(&actors->group_ids, int64_t)) {
            if (darray_init(&actors->group_traits, struct vector *)) {
                if (darray_init(&actors->group_buckets, iproc_group_bucket)) {
                    if (refcount_init(&actors->refcount)) {
                        actors->dim = dim;
                        return actors;
                    }
                    darray_deinit(&actors->group_buckets);
                }
                darray_deinit(&actors->group_traits);
            }
            darray_deinit(&actors->group_ids);
        }
        free(actors);
    }
    return NULL;
}



iproc_actors *
iproc_actors_ref (iproc_actors *actors)
{
    if (actors) {
        refcount_get(&actors->refcount);
    }
    return actors;
}

static void
iproc_actors_release (struct refcount *refcount)
{
    iproc_actors *actors = container_of(refcount, iproc_actors, refcount);
    iproc_actors_free(actors);
}

void
iproc_actors_unref (iproc_actors *actors)
{
    if (!actors)
        return;

    refcount_put(&actors->refcount, iproc_actors_release);
}


ssize_t iproc_actors_size (const iproc_actors *actors)
{
    assert(actors);
    return darray_size(&actors->group_ids);
}

ssize_t iproc_actors_dim (const iproc_actors *actors)
{
    assert(actors);
    return actors->dim;
}

bool iproc_actors_reserve (iproc_actors *actors, ssize_t size)
{
    assert(actors);
    assert(size >= 0);
    
    return darray_reserve(&actors->group_ids, size);
}


ssize_t iproc_actors_add (iproc_actors *actors, const struct vector *traits)
{
    assert(actors);
    assert(vector_size(traits) == iproc_actors_dim(actors));

    ssize_t id = -1;
    
    if (iproc_actors_reserve(actors, iproc_actors_size(actors) + 1)) {
        id = iproc_actors_size(actors);
    
        int64_t group_id = iproc_actors_insert_group(actors, traits);
    
        darray_push_back(&actors->group_ids, &group_id);
    }
        
    return id;
}


struct vector *
iproc_actors_get (iproc_actors *actors,
                  int64_t       actor_id)
{
    assert(actors);
    assert(0 <= actor_id);
    assert(actor_id < iproc_actors_size(actors));

    int64_t g = iproc_actors_group(actors, actor_id);
    struct vector *x = iproc_actors_group_traits(actors, g);
    return x;
}

int64_t
iproc_actors_ngroup (iproc_actors *actors)
{
    assert(actors);
    return darray_size(&actors->group_traits);
}

int64_t
iproc_actors_group (iproc_actors *actors,
                    int64_t       actor_id)
{
    assert(actors);
    assert(0 <= actor_id);
    assert(actor_id < iproc_actors_size(actors));

    struct darray *group_ids = &actors->group_ids;
    int64_t g = darray_index(group_ids, int64_t, actor_id);
    return g;
}


struct vector *
iproc_actors_group_traits (iproc_actors *actors,
                           int64_t       group_id)
{
    assert(actors);
    assert(0 <= group_id);
    assert(group_id < iproc_actors_ngroup(actors));

    struct vector *x = darray_index(&actors->group_traits,
                                   struct vector *,
                                   group_id);
    return x;
}


void
iproc_actors_mul (double        alpha,
                  iproc_trans   trans,
                  iproc_actors *actors,
                  struct vector *x,
                  double        beta,
                  struct vector *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || vector_size(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || vector_size(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || vector_size(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || vector_size(y) == iproc_actors_dim(actors));

    int64_t n = iproc_actors_size(actors);
    int64_t i;
    struct vector *row;
    double dot, entry;

    if (beta == 0) {
        vector_fill(y, 0.0);
    } else if (beta != 1) {
        vector_scale(y, beta);
    }

    /* NOTE: this could be more efficient by using info about actor groups */
    if (trans == IPROC_TRANS_NOTRANS) {
        for (i = 0; i < n; i++) {
            row = iproc_actors_get(actors, i);
            dot = vector_dot(row, x);
            vector_index(y, i) += alpha * dot;
        }
    } else {
        for (i = 0; i < n; i++) {
            row = iproc_actors_get(actors, i);
            entry = vector_index(x, i);
            vector_axpy(alpha * entry, row, y);
        }
    }
}


void
iproc_actors_muls (double         alpha,
                   iproc_trans    trans,
                   iproc_actors  *actors,
                   iproc_svector *x,
                   double         beta,
                   struct vector  *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || vector_size(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_svector_dim(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || vector_size(y) == iproc_actors_dim(actors));

    int64_t n = iproc_actors_size(actors);
    int64_t i;
    struct vector *row;
    double dot, entry;

    if (beta == 0) {
        vector_fill(y, 0.0);
    } else if (beta != 1) {
        vector_scale(y, beta);
    }

    /* NOTE: this could be more efficient by using info about actor groups */
    if (trans == IPROC_TRANS_NOTRANS) {
        for (i = 0; i < n; i++) {
            row = iproc_actors_get(actors, i);
            dot = iproc_vector_sdot(row, x);
            vector_index(y, i) += alpha * dot;
        }
    } else {
        int64_t inz, nnz = iproc_svector_nnz(x);
        for (inz = 0; inz < nnz; inz++) {
            i = iproc_svector_nz(x, inz);
            row = iproc_actors_get(actors, i);
            entry = iproc_svector_nz_get(x, inz);
            vector_axpy(alpha * entry, row, y);
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

    struct vector xcol;
    struct vector ycol;
    
    for (j = 0; j < m; j++) {
        vector_init_matrix_col(&xcol, x, j);
        vector_init_matrix_col(&ycol, y, j);
        iproc_actors_mul(alpha, trans, actors, &xcol, 1.0, &ycol);
    }
}
