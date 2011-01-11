#include <assert.h>
#include <stdlib.h>
#include <strings.h>
#include "blas-private.h"
#include "memory.h"
#include "utils.h"
#include "matrix.h"


static void
iproc_matrix_free (iproc_matrix *matrix)
{
    if (matrix) {
        iproc_free(matrix->data);
        iproc_free(matrix);
    }
}

iproc_matrix *
iproc_matrix_new (int64_t nrow, int64_t ncol)
{
    iproc_matrix *matrix = iproc_malloc(sizeof(*matrix));

    if (!matrix) return NULL;

    matrix->data = iproc_malloc(nrow * ncol * sizeof(matrix->data[0]));
    matrix->nrow = nrow;
    matrix->ncol = ncol;
    matrix->lda = IPROC_MAX(1, nrow);
    iproc_refcount_init(&matrix->refcount);

    if (!(matrix->data)) {
        iproc_matrix_free(matrix);
        matrix = NULL;
    }

    return matrix;
}

iproc_matrix *
iproc_matrix_new_copy (iproc_matrix *matrix)
{
    assert(matrix);
    int64_t m = iproc_matrix_nrow(matrix);
    int64_t n = iproc_matrix_ncol(matrix);
    iproc_matrix *copy = iproc_matrix_new(m, n);
    iproc_matrix_copy(copy, matrix);
    return copy;
}

iproc_matrix *
iproc_matrix_ref (iproc_matrix *matrix)
{
    if (matrix) {
        iproc_refcount_get(&matrix->refcount);
    }
    return matrix;
}

static void
iproc_matrix_release (iproc_refcount *refcount)
{
    iproc_matrix *matrix = container_of(refcount, iproc_matrix, refcount);
    iproc_matrix_free(matrix);
}

void
iproc_matrix_unref (iproc_matrix *matrix)
{
    if (!matrix)
        return;

    iproc_refcount_put(&matrix->refcount, iproc_matrix_release);
}

int64_t
iproc_matrix_nrow (iproc_matrix *matrix)
{
    assert(matrix);
    return matrix->nrow;
}

int64_t
iproc_matrix_ncol (iproc_matrix *matrix)
{
    assert(matrix);
    return matrix->ncol;
}

int64_t
iproc_matrix_lda (iproc_matrix *matrix)
{
    assert(matrix);
    return matrix->lda;
}

void
iproc_matrix_set_all (iproc_matrix *matrix,
                      double       value)
{
    assert(matrix);

    int64_t m = iproc_matrix_nrow(matrix);
    int64_t n = iproc_matrix_ncol(matrix);
    int64_t i, j;

    for (j = 0; j < n; j++) {
        for (i = 0; i < m; i++) {
            iproc_matrix_set(matrix, i, j, value);
        }
    }
}

void
iproc_matrix_copy (iproc_matrix *dst_matrix,
                   iproc_matrix *matrix)
{
    assert(dst_matrix);
    assert(matrix);
    assert(iproc_matrix_nrow(dst_matrix) == iproc_matrix_nrow(matrix));
    assert(iproc_matrix_ncol(dst_matrix) == iproc_matrix_ncol(matrix));

    int64_t m = iproc_matrix_nrow(dst_matrix);
    int64_t n = iproc_matrix_ncol(dst_matrix);
    int64_t i, j;

    if (iproc_matrix_lda(dst_matrix) == m && iproc_matrix_lda(matrix) == m) {
        int64_t mn = m * n;
        int64_t one = 1;

        F77_FUNC(dcopy)(F77_INTP(mn), iproc_matrix_ptr(matrix, 0, 0),
                        F77_INTP(one), iproc_matrix_ptr(dst_matrix, 0, 0),
                        F77_INTP(one));
    } else {
        double value;

        for (j = 0; j < n; j++) {
            for (i = 0; i < m; i++) {
                value = iproc_matrix_get(matrix, i, j);
                iproc_matrix_set(dst_matrix, i, j, value);
            }
        }
    }
}


double
iproc_matrix_get (iproc_matrix *a, int64_t i, int64_t j)
{
    assert(a);
    assert(0 <= i && i < iproc_matrix_nrow(a));
    assert(0 <= j && j < iproc_matrix_ncol(a));

    double *ptr = iproc_matrix_ptr(a, i, j);
    return *ptr;
}

void
iproc_matrix_set (iproc_matrix *a, int64_t i, int64_t j,
          double val)
{
    assert(a);
    assert(0 <= i && i < iproc_matrix_nrow(a));
    assert(0 <= j && j < iproc_matrix_ncol(a));

    double *ptr = iproc_matrix_ptr(a, i, j);
    *ptr = val;
}

double *
iproc_matrix_ptr (iproc_matrix *a, int64_t i, int64_t j)
{
    assert(a);
    assert(0 <= i && i <= iproc_matrix_nrow(a));
    assert(0 <= j && j <= iproc_matrix_ncol(a));

    return a->data + (i + j * (iproc_matrix_lda(a)));
}

void
iproc_matrix_add (iproc_matrix *dst_matrix,
                  iproc_matrix *matrix)
{
    iproc_matrix_acc(dst_matrix, +1.0, matrix);
}

void
iproc_matrix_sub (iproc_matrix *dst_matrix,
                  iproc_matrix *matrix)
{
    iproc_matrix_acc(dst_matrix, -1.0, matrix);
}

void
iproc_matrix_acc (iproc_matrix *dst_matrix,
                  double        scale,
                  iproc_matrix *matrix)
{
    assert(dst_matrix);
    assert(matrix);
    assert(iproc_matrix_nrow(dst_matrix) == iproc_matrix_nrow(matrix));
    assert(iproc_matrix_ncol(dst_matrix) == iproc_matrix_ncol(matrix));

    int64_t m = iproc_matrix_nrow(dst_matrix);
    int64_t n = iproc_matrix_ncol(dst_matrix);
    int64_t one = 1;
    int64_t j;

    if (iproc_matrix_lda(matrix) == m && iproc_matrix_lda(dst_matrix) == m) {
        int64_t mn = m * n;
        F77_FUNC(daxpy)(F77_INTP(mn), &scale,
                        iproc_matrix_ptr(matrix, 0, 0), F77_INTP(one),
                        iproc_matrix_ptr(dst_matrix, 0, 0), F77_INTP(one));
    } else {
        for (j = 0; j < n; j++) {
            F77_FUNC(daxpy)(F77_INTP(m), &scale,
                            iproc_matrix_ptr(matrix, 0, j), F77_INTP(one),
                            iproc_matrix_ptr(dst_matrix, 0, j), F77_INTP(one));
        }
    }

}

void
iproc_matrix_scale (iproc_matrix *matrix,
                    double        scale)
{
    assert(matrix);

    int64_t n = iproc_matrix_ncol(matrix);
    int64_t j;

    for (j = 0; j < n; j++) {
        iproc_vector_view col = iproc_matrix_col(matrix, j);
        iproc_vector_scale(&col.vector, scale);
    }
}

void
iproc_matrix_scale_rows (iproc_matrix *matrix,
                         iproc_vector *scale)
{
    assert(matrix);
    assert(scale);
    assert(iproc_matrix_nrow(matrix) == iproc_vector_dim(scale));

    int64_t n = iproc_matrix_ncol(matrix);
    int64_t j;

    for (j = 0; j < n; j++) {
        iproc_vector_view col = iproc_matrix_col(matrix, j);
        iproc_vector_mul(&col.vector, scale);
    }
}

iproc_vector_view
iproc_matrix_col (iproc_matrix *matrix,
                  int64_t       j)
{
    assert(0 <= j && j < iproc_matrix_ncol(matrix));
    int64_t dim = iproc_matrix_nrow(matrix);
    double *data = iproc_matrix_ptr(matrix, 0, j);
    return iproc_vector_view_array(data, dim);
}

iproc_matrix_view
iproc_matrix_submatrix (iproc_matrix *matrix,
                        int64_t       i,
                        int64_t       j,
                        int64_t       nrow,
                        int64_t       ncol)
{
    int64_t lda = iproc_matrix_lda(matrix);
    double *data = iproc_matrix_ptr(matrix, i, j);
    return iproc_matrix_view_array_with_lda(data, nrow, ncol, lda);
}

iproc_matrix_view
iproc_matrix_cols (iproc_matrix *matrix,
                   int64_t       j,
                   int64_t       n)
{
    assert(matrix);
    assert(j + n <= iproc_matrix_ncol(matrix));

    int64_t m = iproc_matrix_nrow(matrix);
    return iproc_matrix_submatrix(matrix, 0, j, m, n);
}

iproc_matrix_view
iproc_matrix_view_array (double *data,
                         int64_t nrow,
                         int64_t ncol)
{
    assert(data);
    int64_t lda = IPROC_MAX(1, nrow);
    return iproc_matrix_view_array_with_lda(data, nrow, ncol, lda);
}

iproc_matrix_view
iproc_matrix_view_array_with_lda (double *data,
                                  int64_t nrow,
                                  int64_t ncol,
                                  int64_t lda)
{
    assert(lda >= IPROC_MAX(1, nrow));

    iproc_matrix_view view = {{ data, nrow, ncol, lda, {0} }};
    return view;
}

iproc_matrix_view
iproc_matrix_view_vector (iproc_vector *vector,
                          int64_t       nrow,
                          int64_t       ncol)
{
    assert(vector);
    assert(iproc_vector_dim(vector) == nrow * ncol);
    double *data = iproc_vector_ptr(vector, 0);
    return iproc_matrix_view_array(data, nrow, ncol);
}


void
iproc_matrix_mul (double        alpha,
                  iproc_trans   trans,
                  iproc_matrix *matrix,
                  iproc_vector *x,
                  double        beta,
                  iproc_vector *y)
{
    assert(matrix);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_matrix_ncol(matrix));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_matrix_nrow(matrix));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_matrix_nrow(matrix));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_matrix_ncol(matrix));

    char    *ptrans = (trans == IPROC_TRANS_NOTRANS) ? "N" : "T";
    int64_t  m      = iproc_matrix_nrow(matrix);
    int64_t  n      = iproc_matrix_ncol(matrix);
    void    *pa     = iproc_matrix_ptr(matrix, 0, 0);
    int64_t  lda    = iproc_matrix_lda(matrix);
    void    *px     = iproc_vector_ptr(x, 0);
    int64_t  incx   = 1;
    void    *py     = iproc_vector_ptr(y, 0);
    int64_t  incy   = 1;

    F77_FUNC(dgemv)(ptrans, F77_INTP(m),  F77_INTP(n), &alpha, pa, F77_INTP(lda), 
                    px, F77_INTP(incx), &beta, py, F77_INTP(incy));
}

void
iproc_matrix_matmul (double        alpha,
                     iproc_trans   trans,
                     iproc_matrix *matrix,
                     iproc_matrix *x,
                     double        beta,
                     iproc_matrix *y)
{
    assert(matrix);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(x) == iproc_matrix_ncol(matrix));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(y) == iproc_matrix_nrow(matrix));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(x) == iproc_matrix_nrow(matrix));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(y) == iproc_matrix_ncol(matrix));
    assert(iproc_matrix_ncol(x) == iproc_matrix_ncol(y));

    char    *ptransa = (trans == IPROC_TRANS_NOTRANS) ? "N" : "T";
    char    *ptransb = "N";
    int64_t  m       = iproc_matrix_nrow(y);
    int64_t  n       = iproc_matrix_ncol(y);
    int64_t  k       = (trans == IPROC_TRANS_NOTRANS) ? iproc_matrix_ncol(matrix)
                                                      : iproc_matrix_nrow(matrix);
    void    *pa      = iproc_matrix_ptr(matrix, 0, 0);
    int64_t  lda     = iproc_matrix_lda(matrix);
    void    *pb      = iproc_matrix_ptr(x, 0, 0);
    int64_t  ldb     = iproc_matrix_lda(x);
    void    *pc      = iproc_matrix_ptr(y, 0, 0);
    int64_t  ldc     = iproc_matrix_lda(y);

    F77_FUNC(dgemm)(ptransa, ptransb, F77_INTP(m), F77_INTP(n), F77_INTP(k),
                    &alpha, pa, F77_INTP(lda), pb, F77_INTP(ldb),
                    &beta, pc, F77_INTP(ldc));
}
