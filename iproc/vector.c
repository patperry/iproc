
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include <math.h>
#include <iproc/memory.h>
#include <iproc/blas-private.h>
#include <iproc/vector.h>


iproc_vector *
iproc_vector_new (int64_t  dim)
{
    iproc_vector *vector = NULL;
    size_t        size   = dim * sizeof(double);

    if (!(dim >= 0)) return NULL;
    if (!(dim <= F77_INT_MAX)) return NULL;
    if (!(dim <= SIZE_MAX / sizeof(double))) return NULL;

    vector = iproc_malloc(sizeof(iproc_vector));
    vector->pdata = iproc_malloc(size);
    vector->dim = dim;
    vector->refcount = 1;

    return vector;
}


static void
iproc_vector_free (iproc_vector *vector)
{
    if (vector) {
        iproc_free(vector->pdata);
        iproc_free(vector);
    }
}

iproc_vector *
iproc_vector_ref (iproc_vector *vector)
{
    if (vector) {
        vector->refcount = vector->refcount + 1;
    }

    return vector;
}

void
iproc_vector_unref (iproc_vector *vector)
{
    if (!vector)
        return;

    if (vector->refcount == 1) {
        iproc_vector_free(vector);
    } else {
        vector->refcount = vector->refcount - 1;
    }
}

int64_t
iproc_vector_dim (iproc_vector *vector)
{
    assert(vector);
    return vector->dim;
}

void
iproc_vector_set_all (iproc_vector *vector,
                      double        value)
{
    assert(vector);

    int64_t n   = iproc_vector_dim(vector);
    double *ptr = iproc_vector_ptr(vector, 0);
    double *end = iproc_vector_ptr(vector, n);

    while (ptr != end) {
        *ptr++ = value;
    }
}


void
iproc_vector_set_basis (iproc_vector *vector,
                        int64_t       index)
{
    assert(vector);
    assert(iproc_vector_dim(vector) > 0);
    assert(0 <= index);
    assert(index < iproc_vector_dim(vector));

    iproc_vector_set_all(vector, 0);
    iproc_vector_set(vector, index, 1);
}


double
iproc_vector_get (iproc_vector *vector,
                  int64_t       index)
{
    assert(vector);
    assert(0 <= index);
    assert(index < iproc_vector_dim(vector));

    double *ptr = iproc_vector_ptr(vector, index);
    return *ptr;
}


void
iproc_vector_set (iproc_vector *vector,
                  int64_t       index,
                  double        value)
{
    assert(vector);
    assert(0 <= index);
    assert(index < iproc_vector_dim(vector));

    double *ptr = iproc_vector_ptr(vector, index);
    *ptr = value;
}

void
iproc_vector_inc (iproc_vector *vector,
                  int64_t       index,
                  double        value)
{
    assert(vector);
    assert(0 <= index);
    assert(index < iproc_vector_dim(vector));

    double *ptr = iproc_vector_ptr(vector, index);
    *ptr += value;
}

double *
iproc_vector_ptr (iproc_vector *vector,
                  int64_t       index)
{
    assert(vector);

    double *ptr = vector->pdata + index;
    return ptr;
}


iproc_vector_view
iproc_vector_subvector (iproc_vector *vector,
                        int64_t       index,
                        int64_t       dim)
{
    assert(vector);
    assert(0 <= index);
    assert(index <= iproc_vector_dim(vector));
    assert(0 <= dim);
    assert(dim <= iproc_vector_dim(vector) - index);

    double *ptr = iproc_vector_ptr(vector, index);
    iproc_vector_view view = iproc_vector_view_array(ptr, dim);
    return view;
}


iproc_vector_view
iproc_vector_view_array (double  *array,
                         int64_t  dim)
{
    assert(dim >= 0);
    assert(array || dim == 0);

    iproc_vector_view view = {{ array, dim, 0 }};
    return view;
}


void
iproc_vector_copy (iproc_vector *dst_vector,
                   iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    int64_t n    = iproc_vector_dim(dst_vector);
    double *px   = iproc_vector_ptr(vector, 0);
    int64_t incx = 1;
    double *py   = iproc_vector_ptr(dst_vector, 0);
    int64_t incy = 1;

    F77_FUNC(dcopy)(F77_INTP(n), px, F77_INTP(incx), py, F77_INTP(incy));
}


void
iproc_vector_swap (iproc_vector *vector1,
                   iproc_vector *vector2)
{
    assert(vector1);
    assert(vector2);
    assert(iproc_vector_dim(vector1) == iproc_vector_dim(vector2));

    int64_t n    = iproc_vector_dim(vector1);
    double *px   = iproc_vector_ptr(vector1, 0);
    int64_t incx = 1;
    double *py   = iproc_vector_ptr(vector2, 0);
    int64_t incy = 1;

    F77_FUNC(dswap)(F77_INTP(n), px, F77_INTP(incx), py, F77_INTP(incy));
}


void
iproc_vector_swap_elems (iproc_vector *vector,
                         int64_t       index1,
                         int64_t       index2)
{
    assert(vector);
    assert(0 <= index1);
    assert(index1 < iproc_vector_dim(vector));
    assert(0 <= index2);
    assert(index2 < iproc_vector_dim(vector));

    double e1, e2;

    if (index1 != index2) {
        e1 = iproc_vector_get(vector, index1);
        e2 = iproc_vector_get(vector, index2);
        iproc_vector_set(vector, index2, e1);
        iproc_vector_set(vector, index1, e2);
    }
}


void
iproc_vector_reverse (iproc_vector *vector)
{
    assert(vector);

    int64_t n = iproc_vector_dim(vector);
    int64_t i;

    for (i = 0; i < n / 2; i++) {
        iproc_vector_swap_elems(vector, i, n - 1 - i);
    }
}


void
iproc_vector_scale (iproc_vector *vector,
                    double        scale)
{
    assert(vector);

    int64_t n     = iproc_vector_dim(vector);
    double  alpha = scale;
    void   *px    = iproc_vector_ptr(vector, 0);
    int64_t incx  = 1;

    F77_FUNC(dscal)(F77_INTP(n), &alpha, px, F77_INTP(incx));
}


void
iproc_vector_shift (iproc_vector *vector,
                    double        shift)
{
    assert(vector);

    int64_t n   = iproc_vector_dim(vector);
    double *ptr = iproc_vector_ptr(vector, 0);
    double *end = iproc_vector_ptr(vector, n);

    while (ptr != end) {
        *ptr++ += shift;
    }
}


void
iproc_vector_add (iproc_vector *dst_vector,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    iproc_vector_acc(dst_vector, 1, vector);
}


void
iproc_vector_sub (iproc_vector *dst_vector,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    iproc_vector_acc(dst_vector, -1, vector);
}


void
iproc_vector_mul (iproc_vector *dst_vector,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    int64_t n     = iproc_vector_dim(dst_vector);
    int64_t k     = 0;
    double *px    = iproc_vector_ptr(vector, 0);
    int64_t incx  = 1;
    double *py    = iproc_vector_ptr(dst_vector, 0);
    int64_t incy  = 1;

    F77_FUNC(dtbmv)("U", "N", "N", F77_INTP(n), F77_INTP(k),
                    px, F77_INTP(incx), py, F77_INTP(incy));
}


void
iproc_vector_div (iproc_vector *dst_vector,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    int64_t n     = iproc_vector_dim(dst_vector);
    int64_t k     = 0;
    double *px    = iproc_vector_ptr(vector, 0);
    int64_t incx  = 1;
    double *py    = iproc_vector_ptr(dst_vector, 0);
    int64_t incy  = 1;

    F77_FUNC(dtbsv)("U", "N", "N", F77_INTP(n), F77_INTP(k),
                    px, F77_INTP(incx), py, F77_INTP(incy));
}


void
iproc_vector_acc (iproc_vector *dst_vector,
                  double        scale,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    int64_t n     = iproc_vector_dim(dst_vector);
    double  alpha = scale;
    double *px    = iproc_vector_ptr(vector, 0);
    int64_t incx  = 1;
    double *py    = iproc_vector_ptr(dst_vector, 0);
    int64_t incy  = 1;

    F77_FUNC(daxpy)(F77_INTP(n), &alpha, px, F77_INTP(incx),
                    py, F77_INTP(incy));
}


double
iproc_vector_dot (iproc_vector *vector1,
                  iproc_vector *vector2)
{
    assert(vector1);
    assert(vector2);
    assert(iproc_vector_dim(vector1) == iproc_vector_dim(vector2));

    int64_t n    = iproc_vector_dim(vector1);
    double *px   = iproc_vector_ptr(vector1, 0);
    int64_t incx = 1;
    double *py   = iproc_vector_ptr(vector2, 0);
    int64_t incy = 1;

    double dot = F77_FUNC(ddot)(F77_INTP(n), px, F77_INTP(incx), py,
                                F77_INTP(incy));
    return dot;
}


double
iproc_vector_norm (iproc_vector *vector)
{
    assert(vector);

    int64_t n    = iproc_vector_dim(vector);
    void   *px   = iproc_vector_ptr(vector, 0);
    int64_t incx = 1;

    double norm = F77_FUNC(dnrm2)(F77_INTP(n), px, F77_INTP(incx));
    return norm;
}


double
iproc_vector_sum_abs (iproc_vector *vector)
{
    assert(vector);

    int64_t n    = iproc_vector_dim(vector);
    void   *px   = iproc_vector_ptr(vector, 0);
    int64_t incx = 1;

    double sum_abs = F77_FUNC(dasum)(F77_INTP(n), px, F77_INTP(incx));
    return sum_abs;
}


double
iproc_vector_max_abs (iproc_vector *vector)
{
    assert(vector);

    double max_abs = 0;
    double e;
    int64_t i;

    if (iproc_vector_dim(vector) > 0) {
        i = iproc_vector_max_abs_index(vector);
        e = iproc_vector_get(vector, i);
        max_abs = fabs(e);
    }

    return max_abs;
}


int64_t
iproc_vector_max_abs_index (iproc_vector *vector)
{
    assert(vector);
    if (!(iproc_vector_dim(vector) > 0)) return -1;

    int64_t n    = iproc_vector_dim(vector);
    void   *px   = iproc_vector_ptr(vector, 0);
    int64_t incx = 1;

    int64_t index1 = (int64_t)(F77_FUNC(idamax)(F77_INTP(n), px,
                                                F77_INTP(incx)));
    int64_t index  = index1 - 1;
    return index;
}



