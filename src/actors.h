#ifndef _IPROC_ACTORS_H
#define _IPROC_ACTORS_H


#include <stddef.h>
#include <stdint.h>

#include "darray.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"

/*
struct actor {
    struct cohort *cohort;
};

struct cohort {
    struct vector traits;
    struct intset actors;
};

struct _iproc_actors {
    struct darray  actors;
    struct hashset cohorts;
    iproc_refcount refcount;
};
*/



typedef struct _iproc_actors       iproc_actors;
typedef struct _iproc_group        iproc_group;
typedef struct _iproc_group_bucket iproc_group_bucket;

struct _iproc_group {
    struct vector *traits;     /* traits must be the first member */
    int64_t       id;
};

struct _iproc_group_bucket {
    size_t        traits_hash; /* traits_hash must be the first member */
    struct darray groups;
};

struct _iproc_actors {
    struct darray   group_ids;
    struct darray   group_traits;
    struct darray   group_buckets;
    struct refcount refcount;
};

/* makes a copy of traits0 */
iproc_actors * iproc_actors_new          (int64_t       size,
                                          struct vector *traits0);
iproc_actors * iproc_actors_ref          (iproc_actors *actors);
void           iproc_actors_unref        (iproc_actors *actors);


int64_t        iproc_actors_size         (iproc_actors *actors);
int64_t        iproc_actors_dim          (iproc_actors *actors);

/* makes a copy of traits */
void           iproc_actors_set          (iproc_actors *actors,
                                          int64_t       actor_id,
                                          struct vector *traits);
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
