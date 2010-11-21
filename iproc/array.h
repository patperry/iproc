#ifndef _IPROC_ARRAY_H
#define _IPROC_ARRAY_H

#include <stddef.h>
#include <stdint.h>

typedef struct _iproc_array iproc_array;

struct _iproc_array {
    int64_t n;
    int64_t n_max;
    void   *data;
    size_t  elem_size;
};

iproc_array * iproc_array_new       (size_t       elem_size);
void          iproc_array_free      (iproc_array *array);

size_t        iproc_array_elem_size (iproc_array *array);
void          iproc_array_set_size  (iproc_array *array,
                                     int64_t      n);

int64_t       iproc_array_size      (iproc_array *array);
void          iproc_array_set       (iproc_array *array,
                                     int64_t      i,
                                     void        *pe);
void          iproc_array_append    (iproc_array *array,
                                     void        *pe);
void          iproc_array_prepend   (iproc_array *array,
                                     void        *pe);
void          iproc_array_insert    (iproc_array *array,
                                     int64_t      i,
                                     void        *pe);
void          iproc_array_remove    (iproc_array *array,
                                     int64_t      i);

int64_t       iproc_array_lfind     (iproc_array *array,
                                     void        *value,
                                     int        (*compare) (void *, void *));
int64_t       iproc_array_bsearch   (iproc_array *array,
                                     void        *value,
                                     int        (*compare) (void *, void *));

#define       iproc_array_index(a,t,i) (((t *)((a)->data))[(i)])



#endif /* _IPROC_ARRAY_H */
