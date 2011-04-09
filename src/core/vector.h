#ifndef _VECTOR_H
#define _VECTOR_H

#include <stddef.h>
#include "array.h"

typedef struct _iproc_vector_view iproc_vector_view;

struct vector {
    struct array array;
};

struct _iproc_vector_view {
    struct vector vector;
};

/* create, destroy */
struct vector * vector_init       (struct vector *v, ssize_t n);
struct vector * vector_init_view  (struct vector *v, double *ptr, ssize_t n);
struct vector * vector_init_slice (struct vector *v,
                                   const struct vector *parent,
                                   ssize_t i, ssize_t n);
struct vector * vector_init_copy (struct vector *v, const struct vector *src);
void            vector_deinit    (struct vector *v);

struct vector * vector_new      (ssize_t n);
struct vector * vector_new_copy (const struct vector *v);
void            vector_free     (struct vector *v);


/* index */
#define vector_index(v,i) array_index(&(v)->array, double, i)


/* informative */
static inline ssize_t vector_size (const struct vector *v);


/* operations */


/* iteration */
static inline double * vector_begin (const struct vector *v);
static inline double * vector_ptr   (const struct vector *v, ssize_t i);
static inline double * vector_end   (const struct vector *v);


void              vector_fill       (struct vector *vector,
                                              double        value);
void              vector_set_basis     (struct vector *vector,
                                              ssize_t       index);

iproc_vector_view vector_slice     (struct vector *vector,
                                              ssize_t       index,
                                              ssize_t       dim);
iproc_vector_view iproc_vector_view_array    (double       *array,
                                              ssize_t       dim);
void              vector_copy          (struct vector *dst_vector,
                                              const struct vector *vector);
void              vector_swap          (struct vector *vector1,
                                              struct vector *vector2);
void              vector_swap_elems    (struct vector *vector,
                                              ssize_t       i,
                                              ssize_t       j);
void              vector_reverse       (struct vector *vector);
void              vector_scale         (struct vector *vector,
                                              double        scale);
void              vector_shift         (struct vector *vector,
                                              double        shift);
void              vector_add           (struct vector *dst_vector,
                                              struct vector *vector);
void              vector_sub           (struct vector *dst_vector,
                                              struct vector *vector);
void              vector_mul           (struct vector *dst_vector,
                                              struct vector *vector);
void              vector_div           (struct vector *dst_vector,
                                              struct vector *vector);
void              vector_acc           (struct vector *dst_vector,
                                              double        scale,
                                              struct vector *vector);
double            vector_dot           (struct vector *vector1,
                                              struct vector *vector2);
double            vector_norm          (struct vector *vector);

double            vector_sum_abs       (struct vector *vector);
double            vector_max_abs       (struct vector *vector);
ssize_t           vector_max_abs_index (struct vector *vector);
double            vector_max           (struct vector *vector);
ssize_t           vector_max_index     (struct vector *vector);
double            vector_log_sum_exp   (struct vector *vector);

void              vector_exp           (struct vector *vector);

void              vector_printf        (struct vector *vector);

size_t            vector_hash          (struct vector *vector);
int               vector_identical     (struct vector *vector1,
                                              struct vector *vector2);
int               vector_compare       (const void *x1,  const void *x2);
int               vector_ptr_compare   (const void *px1, const void *px2);


/* inline function definitions */
ssize_t vector_size   (const struct vector *v) { return array_size(&v->array); }
double * vector_begin (const struct vector *v) { return vector_ptr(v, 0); }
double * vector_ptr   (const struct vector *v, ssize_t i) { return &vector_index(v, i); }
double * vector_end   (const struct vector *v) { return vector_ptr(v, vector_size(v)); }






#endif /* _VECTOR_H */
