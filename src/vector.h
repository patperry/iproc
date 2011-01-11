#ifndef _IPROC_VECTOR_H
#define _IPROC_VECTOR_H

#include <stdint.h>
#include "refcount.h"

typedef struct _iproc_vector      iproc_vector;
typedef struct _iproc_vector_view iproc_vector_view;

struct _iproc_vector {
    double        *pdata;
    int64_t        dim;
    iproc_refcount refcount;
};

struct _iproc_vector_view {
    iproc_vector vector;
};


iproc_vector *    iproc_vector_new           (int64_t       dim);
iproc_vector *    iproc_vector_new_copy      (iproc_vector *vector);
iproc_vector *    iproc_vector_ref           (iproc_vector *vector);
void              iproc_vector_unref         (iproc_vector *vector);
int64_t           iproc_vector_dim           (iproc_vector *vector);
void              iproc_vector_set_all       (iproc_vector *vector,
                                              double        value);
void              iproc_vector_set_basis     (iproc_vector *vector,
                                              int64_t       index);
double            iproc_vector_get           (iproc_vector *vector,
                                              int64_t       index);
void              iproc_vector_set           (iproc_vector *vector,
                                              int64_t       index,
                                              double        value);
void              iproc_vector_inc           (iproc_vector *vector,
                                              int64_t       index,
                                              double        value);
double *          iproc_vector_ptr           (iproc_vector *vector,
                                              int64_t       index);
iproc_vector_view iproc_vector_subvector     (iproc_vector *vector,
                                              int64_t       index,
                                              int64_t       dim);
iproc_vector_view iproc_vector_view_array    (double       *array,
                                              int64_t       dim);
void              iproc_vector_copy          (iproc_vector *dst_vector,
                                              iproc_vector *vector);
void              iproc_vector_swap          (iproc_vector *vector1,
                                              iproc_vector *vector2);
void              iproc_vector_swap_elems    (iproc_vector *vector,
                                              int64_t       index1,
                                              int64_t       index2);
void              iproc_vector_reverse       (iproc_vector *vector);
void              iproc_vector_scale         (iproc_vector *vector,
                                              double        scale);
void              iproc_vector_shift         (iproc_vector *vector,
                                              double        shift);
void              iproc_vector_add           (iproc_vector *dst_vector,
                                              iproc_vector *vector);
void              iproc_vector_sub           (iproc_vector *dst_vector,
                                              iproc_vector *vector);
void              iproc_vector_mul           (iproc_vector *dst_vector,
                                              iproc_vector *vector);
void              iproc_vector_div           (iproc_vector *dst_vector,
                                              iproc_vector *vector);
void              iproc_vector_acc           (iproc_vector *dst_vector,
                                              double        scale,
                                              iproc_vector *vector);
double            iproc_vector_dot           (iproc_vector *vector1,
                                              iproc_vector *vector2);
double            iproc_vector_norm          (iproc_vector *vector);

double            iproc_vector_sum_abs       (iproc_vector *vector);
double            iproc_vector_max_abs       (iproc_vector *vector);
int64_t           iproc_vector_max_abs_index (iproc_vector *vector);
double            iproc_vector_max           (iproc_vector *vector);
int64_t           iproc_vector_max_index     (iproc_vector *vector);
double            iproc_vector_log_sum_exp   (iproc_vector *vector);

void              iproc_vector_exp           (iproc_vector *vector);


#endif /* _IPROC_VECTOR_H */
