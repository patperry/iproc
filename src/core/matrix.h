#ifndef _IPROC_MATRIX_H
#define _IPROC_MATRIX_H

#include <stdint.h>
#include "refcount.h"
#include "vector.h"

typedef struct _iproc_matrix      iproc_matrix;
typedef struct _iproc_matrix_view iproc_matrix_view;

struct _iproc_matrix {
    double        *data;
    int64_t        nrow;
    int64_t        ncol;
    int64_t        lda;
    iproc_refcount refcount;
};

struct _iproc_matrix_view {
    iproc_matrix matrix;
};

typedef enum _iproc_trans {
    IPROC_TRANS_NOTRANS,
    IPROC_TRANS_TRANS,
    IPROC_TRANS_CONJTRANS
} iproc_trans;


iproc_matrix *    iproc_matrix_new                 (int64_t       nrow,
                                                    int64_t       ncol);
iproc_matrix *    iproc_matrix_new_copy            (iproc_matrix *matrix);
iproc_matrix *    iproc_matrix_ref                 (iproc_matrix *matrix);
void              iproc_matrix_unref               (iproc_matrix *matrix);
void              iproc_matrix_copy                (iproc_matrix *dst,
                                                    iproc_matrix *src);
void              iproc_matrix_set_all             (iproc_matrix *matrix,
                                                    double        value);

int64_t           iproc_matrix_nrow                (iproc_matrix *a);
int64_t           iproc_matrix_ncol                (iproc_matrix *a);
int64_t           iproc_matrix_lda                 (iproc_matrix *a);

double            iproc_matrix_get                 (iproc_matrix *a,
                                                    int64_t       i,
                                                    int64_t       j);
void              iproc_matrix_set                 (iproc_matrix *a,
                                                    int64_t       i,
                                                    int64_t       j,
                                                    double        val);
double *          iproc_matrix_ptr                 (iproc_matrix *a,
                                                    int64_t       i,
                                                    int64_t       j);

void              iproc_matrix_add                 (iproc_matrix *dst_matrix,
                                                    iproc_matrix *matrix);
void              iproc_matrix_sub                 (iproc_matrix *dst_matrix,
                                                    iproc_matrix *matrix);
void              iproc_matrix_acc                 (iproc_matrix *dst_matrix,
                                                    double        scale,
                                                    iproc_matrix *matrix);
void              iproc_matrix_scale               (iproc_matrix *matrix,
                                                    double        scale);
void              iproc_matrix_scale_rows          (iproc_matrix *matrix,
                                                    iproc_vector *scale);

iproc_vector_view iproc_matrix_col                 (iproc_matrix *a,
                                                    int64_t       j);
iproc_matrix_view iproc_matrix_cols                (iproc_matrix *matrix,
                                                    int64_t       j,
                                                    int64_t       n);
iproc_matrix_view iproc_matrix_submatrix           (iproc_matrix *a,
                                                    int64_t i,
                                                    int64_t j,
                                                    int64_t nrow,
                                                    int64_t ncol);
iproc_matrix_view iproc_matrix_view_array          (double       *data,
                                                    int64_t       nrow,
                                                    int64_t       ncol);
iproc_matrix_view iproc_matrix_view_array_with_lda (double       *data,
                                                    int64_t       nrow,
                                                    int64_t       ncol,
                                                    int64_t       lda);
iproc_matrix_view iproc_matrix_view_vector         (iproc_vector *vector,
                                                    int64_t       nrow,
                                                    int64_t       ncol);

void              iproc_matrix_mul                 (double        alpha,
                                                    iproc_trans   trans,
                                                    iproc_matrix *matrix,
                                                    iproc_vector *x,
                                                    double        beta,
                                                    iproc_vector *y);
void              iproc_matrix_matmul              (double        alpha,
                                                    iproc_trans   trans,
                                                    iproc_matrix *matrix,
                                                    iproc_matrix *x,
                                                    double        beta,
                                                    iproc_matrix *y);

#endif /* _IPROC_MATRIX_H */
