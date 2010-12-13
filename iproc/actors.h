#ifndef _IPROC_ACTORS_H
#define _IPROC_ACTORS_H

#include <iproc/array.h>
#include <iproc/vector.h>


#define IPROC_ACTORS_DEFCLASS   0

typedef struct _iproc_actors iproc_actors;

struct _iproc_actors {
    iproc_array *vector;
    iproc_array *class;
    int          refcount;
};

iproc_actors * iproc_actors_new          (int64_t       size,
                                          iproc_vector *defvector);
iproc_actors * iproc_actors_ref          (iproc_actors *actors);
void           iproc_actors_unref        (iproc_actors *actors);

int64_t        iproc_actors_append_class (iproc_actors *actors,
                                          iproc_vector *vector);
iproc_vector * iproc_actors_class_vector (iproc_actors *actors,
                                          int64_t       c);
int64_t        iproc_actors_nclass       (iproc_actors *actors);
int64_t        iproc_actors_size         (iproc_actors *actors);
int64_t        iproc_actors_dim          (iproc_actors *actors);


void           iproc_actors_set          (iproc_actors *actors,
                                          int64_t       i,
                                          int64_t       c);
int64_t        iproc_actors_class        (iproc_actors *actors,
                                          int64_t       i);
iproc_vector * iproc_actors_vector       (iproc_actors *actors,
                                          int64_t       i);



#endif /* _IPROC_ACTORS_H */
