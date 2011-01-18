#ifndef _IPROC_ACTORS_H
#define _IPROC_ACTORS_H

#include "array.h"
#include "matrix.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"


#define IPROC_ACTORS_DEFGROUP   0

typedef struct _iproc_actors iproc_actors;

struct _iproc_actors {
    iproc_array   *group_traits;
    iproc_array   *groups;
    iproc_refcount refcount;
};

/* Makes a copy of defvector */
iproc_actors * iproc_actors_new          (int64_t       size,
                                          iproc_vector *traits0);
iproc_actors * iproc_actors_ref          (iproc_actors *actors);
void           iproc_actors_unref        (iproc_actors *actors);

/* Makes a copy of vector */
int64_t        iproc_actors_append_group (iproc_actors *actors,
                                          iproc_vector *traits);
iproc_vector * iproc_actors_group_traits (iproc_actors *actors,
                                          int64_t       g);
int64_t        iproc_actors_ngroup       (iproc_actors *actors);
int64_t        iproc_actors_size         (iproc_actors *actors);
int64_t        iproc_actors_dim          (iproc_actors *actors);

void           iproc_actors_set          (iproc_actors *actors,
                                          int64_t       i,
                                          int64_t       g);
int64_t        iproc_actors_group        (iproc_actors *actors,
                                          int64_t       i);
iproc_vector * iproc_actors_traits       (iproc_actors *actors,
                                          int64_t       i);

void           iproc_actors_mul          (double        alpha,
                                          iproc_trans   trans,
                                          iproc_actors *actors,
                                          iproc_vector *x,
                                          double        beta,
                                          iproc_vector *y);
void           iproc_actors_muls         (double         alpha,
                                          iproc_trans    trans,
                                          iproc_actors  *actors,
                                          iproc_svector *x,
                                          double         beta,
                                          iproc_vector  *y);

void           iproc_actors_matmul       (double        alpha,
                                          iproc_trans   trans,
                                          iproc_actors *actors,
                                          iproc_matrix *x,
                                          double        beta,
                                          iproc_matrix *y);


#endif /* _IPROC_ACTORS_H */
