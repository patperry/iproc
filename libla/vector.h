#ifndef _LA_VECTOR_H
#define _LA_VECTOR_H

#include <libla/index.h>

typedef struct _LAVector     LAVector;
typedef struct _LAVectorView LAVectorView;

struct _LAVector {
    double *pdata;
    la_size dim;
};

struct _LAVectorView {
    LAVector vector;
};


LAVector *   la_vector_new           (la_size   dim);
void         la_vector_free          (LAVector *vector);
la_size      la_vector_dim           (LAVector *vector);
void         la_vector_set_all       (LAVector *vector,
                                      double    value);
double       la_vector_get           (LAVector *vector,
                                      la_index  index);
void         la_vector_set           (LAVector *vector,
                                      la_index  index,
                                      double    value);
double *     la_vector_ptr           (LAVector *vector,
                                      la_index  index);
LAVectorView la_vector_subvector     (LAVector *vector,
                                      la_index index,
                                      la_size  dim);
LAVectorView la_vector_view_array    (double  *array,
                                      la_size  dim);
void         la_vector_copy          (LAVector *dst_vector,
                                      LAVector *vector);
void         la_vector_swap          (LAVector *vector1,
                                      LAVector *vector2);
void         la_vector_swap_elems    (LAVector *vector,
                                      la_index  index1,
                                      la_index  index2);
void         la_vector_reverse       (LAVector *vector);
void         la_vector_add           (LAVector *dst_vector,
                                      LAVector *vector);
void         la_vector_sub           (LAVector *dst_vector,
                                      LAVector *vector);
void         la_vector_mul           (LAVector *dst_vector,
                                      LAVector *vector);
void         la_vector_div           (LAVector *dst_vector,
                                      LAVector *vector);
void         la_vector_acc           (LAVector *dst_vector,
                                      double    scale,
                                      LAVector *vector);
void         la_vector_scale         (LAVector *vector,
                                      double    scale);
void         la_vector_shift         (LAVector *vector,
                                      double    shift);
double       la_vector_dot           (LAVector *vector1,
                                      LAVector *vector2);
double       la_vector_norm          (LAVector *vector);
double       la_vector_sum_abs       (LAVector *vector);
double       la_vector_max_abs       (LAVector *vector);
la_index     la_vector_max_abs_index (LAVector *vector);

#endif /* _LA_VECTOR_H */
