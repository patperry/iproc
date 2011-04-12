#ifndef _IPROC_ACTORS_H
#define _IPROC_ACTORS_H


#include <stddef.h>

#include "darray.h"
#include "hashset.h"
#include "refcount.h"

#include "matrix.h"
#include "svector.h"
#include "vector.h"

typedef struct _iproc_actors iproc_actors;

/*
struct actor {
    struct cohort *cohort;
};

struct cohort {
    struct vector traits;
    struct intset actors;
};

struct _iproc_actors {
    struct hashset  cohorts;
    struct darray   actors;
    ssize_t         dim; 
    struct refcount refcount;
};
*/


typedef struct _iproc_group iproc_group;
struct _iproc_group {
    struct vector *traits;     /* traits must be the first member */
    int64_t       id;
};

typedef struct _iproc_group_bucket iproc_group_bucket;
struct _iproc_group_bucket {
    size_t        traits_hash; /* traits_hash must be the first member */
    struct darray groups;
};

struct _iproc_actors {
    ssize_t         dim;
    struct darray   group_ids;
    struct darray   group_traits;
    struct darray   group_buckets;
    struct refcount refcount;
};

/* makes a copy of traits0 */
iproc_actors * iproc_actors_new   (ssize_t dim);
iproc_actors * iproc_actors_ref   (iproc_actors *actors);
void           iproc_actors_unref (iproc_actors *actors);


ssize_t        iproc_actors_size  (const iproc_actors *actors);
ssize_t        iproc_actors_dim   (const iproc_actors *actors);

/* makes a copy of traits */
ssize_t         iproc_actors_add  (iproc_actors *actors,
                                   const struct vector *traits);
struct vector * iproc_actors_get          (iproc_actors *actors,
                                          int64_t       actor_id);

int64_t        iproc_actors_ngroup       (iproc_actors *actors);
int64_t        iproc_actors_group        (iproc_actors *actors,
                                          int64_t       actor_id);
struct vector * iproc_actors_group_traits (iproc_actors *actors,
                                          int64_t       group_id);


void           iproc_actors_mul          (double        alpha,
                                          iproc_trans   trans,
                                          iproc_actors *actors,
                                          struct vector *x,
                                          double        beta,
                                          struct vector *y);
void           iproc_actors_muls         (double         alpha,
                                          iproc_trans    trans,
                                          iproc_actors  *actors,
                                          iproc_svector *x,
                                          double         beta,
                                          struct vector  *y);

void           iproc_actors_matmul       (double        alpha,
                                          iproc_trans   trans,
                                          iproc_actors *actors,
                                          iproc_matrix *x,
                                          double        beta,
                                          iproc_matrix *y);




#endif /* _IPROC_ACTORS_H */
