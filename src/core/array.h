#ifndef _IPROC_ARRAY_H
#define _IPROC_ARRAY_H

#include <stdbool.h>
#include <sys/types.h>
#include "refcount.h"

typedef struct _iproc_array iproc_array;

struct _iproc_array {
    void          *data;
    ssize_t        n;
    ssize_t        n_max;
    size_t         elem_size;
    iproc_refcount refcount;
};

iproc_array * iproc_array_new       (size_t       elem_size);
iproc_array * iproc_array_new_copy  (iproc_array *array);
void          iproc_array_ref       (iproc_array *array);
void          iproc_array_unref     (iproc_array *array);


size_t        iproc_array_elem_size (iproc_array *array);
void          iproc_array_set_size  (iproc_array *array,
                                     ssize_t      n);

bool          iproc_array_empty     (iproc_array *array);
ssize_t       iproc_array_size      (iproc_array *array);

/* Grow the array if necessary, padding it with zeroes */
void          iproc_array_set_size  (iproc_array *array,
                                     ssize_t      n);

void          iproc_array_set       (iproc_array *array,
                                     ssize_t      i,
                                     void        *pe);
void          iproc_array_append    (iproc_array *array,
                                     void        *pe);
void          iproc_array_prepend   (iproc_array *array,
                                     void        *pe);
void          iproc_array_insert    (iproc_array *array,
                                     ssize_t      i,
                                     void        *pe);
void          iproc_array_remove    (iproc_array *array,
                                     ssize_t      i);

ssize_t       iproc_array_lfind     (iproc_array *array,
                                     void        *value,
                                     int        (*compare) (void *, void *));
ssize_t       iproc_array_bsearch   (iproc_array *array,
                                     void        *value,
                                     int        (*compare) (void *, void *));


#define       iproc_array_index(a,t,i) (((t *)((a)->data))[(i)])



#endif /* _IPROC_ARRAY_H */
