
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>
#include <glib.h>
#include <libla/blas-private.h>
#include <libla/vector.h>


LAVector *
la_vector_new (la_size dim)
{
    LAVector *vector = NULL;
    gsize    size    = dim * sizeof(double);

    g_return_val_if_fail(dim >= 0, NULL);
    g_return_val_if_fail(dim <= G_MAXSIZE / sizeof(double), NULL);

    vector = g_malloc(sizeof(LAVector));
    vector->pdata = g_malloc(size);
    vector->dim = dim;

    return vector;
}


void
la_vector_free (LAVector *vector)
{
    if (vector) {
        g_free(vector->pdata);
        g_free(vector);
    }
}


la_size
la_vector_dim (LAVector *vector)
{
    g_assert(vector);
    return vector->dim;
}


void
la_vector_set_all (LAVector *vector,
                   double    value)
{
    g_assert(vector);

    la_size n   = la_vector_dim(vector);
    double *ptr = la_vector_ptr(vector, 0);
    double *end = la_vector_ptr(vector, n);

    while (ptr != end) {
        *ptr++ = value;
    }
}


double
la_vector_get (LAVector *vector,
               la_index  index)
{
    g_assert(vector);
    g_assert(0 <= index);
    g_assert(index < la_vector_dim(vector));

    double *ptr = la_vector_ptr(vector, index);
    return *ptr;
}


void
la_vector_set (LAVector *vector,
               la_index  index,
               double    value)
{
    g_assert(vector);
    g_assert(0 <= index);
    g_assert(index < la_vector_dim(vector));

    double *ptr = la_vector_ptr(vector, index);
    *ptr = value;
}


double *
la_vector_ptr (LAVector *vector,
               la_index  index)
{
    g_assert(vector);

    double *ptr = vector->pdata + index;
    return ptr;
}


LAVectorView
la_vector_subvector (LAVector *vector,
                     la_index  index,
                     la_size   dim)
{
    g_assert(vector);
    g_assert(0 <= index);
    g_assert(index <= la_vector_dim(vector));
    g_assert(0 <= dim);
    g_assert(dim <= la_vector_dim(vector) - index);

    double *ptr       = la_vector_ptr(vector, index);
    LAVectorView view = la_vector_view_array(ptr, dim);
    return view;
}


LAVectorView
la_vector_view_array (double  *array,
                      la_size  dim)
{
    g_assert(dim >= 0);
    g_assert(array || dim == 0);

    LAVectorView view = {{ array, dim }};
    return view;
}


void
la_vector_copy (LAVector *dst_vector,
LAVector *vector)
{
    g_assert(dst_vector);
    g_assert(vector);
    g_assert(la_vector_dim(dst_vector) == la_vector_dim(vector));

    la_size n    = la_vector_dim(dst_vector);
    double *px   = la_vector_ptr(vector, 0);
    la_size incx = 1;
    double *py   = la_vector_ptr(dst_vector, 0);
    la_size incy = 1;

    F77_FUNC(dcopy)(&n, px, &incx, py, &incy);
}


void
la_vector_swap (LAVector *vector1,
LAVector *vector2)
{
    g_assert(vector1);
    g_assert(vector2);
    g_assert(la_vector_dim(vector1) == la_vector_dim(vector2));

    la_size n    = la_vector_dim(vector1);
    double *px   = la_vector_ptr(vector1, 0);
    la_size incx = 1;
    double *py   = la_vector_ptr(vector2, 0);
    la_size incy = 1;

    F77_FUNC(dswap)(&n, px, &incx, py, &incy);
}


void
la_vector_swap_elems (LAVector *vector,
la_index  index1,
la_index  index2)
{
    g_assert(vector);
    g_assert(0 <= index1);
    g_assert(index1 < la_vector_dim(vector));
    g_assert(0 <= index2);
    g_assert(index2 < la_vector_dim(vector));

    double e1, e2;

    if (index1 != index2) {
        e1 = la_vector_get(vector, index1);
        e2 = la_vector_get(vector, index2);
        la_vector_set_real(vector, index2, e1);
        la_vector_set_real(vector, index1, e2);
    }
}


void
la_vector_reverse (LAVector *vector)
{
    g_assert(vector);

    la_index n = la_vector_dim(vector);
    la_index i;

    for (i = 0; i < n / 2; i++) {
        la_vector_swap_elems(vector, i, n - 1 - i);
    }
}


void
la_vector_scale (LAVector *vector,
double    scale)
{
    g_assert(vector);

    la_size n     = la_vector_dim(vector);
    double  alpha = scale;
    void   *px    = la_vector_ptr(vector, 0);
    la_size incx  = 1;

    F77_FUNC(dscal)(&n, &alpha, px, &incx);
}


void
la_vector_shift (LAVector *vector,
double    shift)
{
    g_assert(vector);

    la_size  n    = la_vector_dim(vector);
    double *ptr = la_vector_ptr(vector, 0);
    double *end = la_vector_ptr(vector, n);

    while (ptr != end) {
        *ptr += shift;
    }
}


void
la_vector_add (LAVector *dst_vector,
LAVector *vector)
{
    g_assert(dst_vector);
    g_assert(vector);
    g_assert(la_vector_dim(dst_vector) == la_vector_dim(vector));

    la_vector_acc(dst_vector, 1, vector);
}


void
la_vector_sub (LAVector *dst_vector,
LAVector *vector)
{
    g_assert(dst_vector);
    g_assert(vector);
    g_assert(la_vector_dim(dst_vector) == la_vector_dim(vector));

    la_vector_acc(dst_vector, -1, vector);
}


void
la_vector_mul (LAVector *dst_vector,
LAVector *vector)
{
    g_assert(dst_vector);
    g_assert(vector);
    g_assert(la_vector_dim(dst_vector) == la_vector_dim(vector));

    la_size n     = la_vector_dim(dst_vector);
    la_size k     = 0;
    double *px    = la_vector_ptr(vector, 0);
    la_size incx  = 1;
    double *py    = la_vector_ptr(dst_vector, 0);
    la_size incy  = 1;

    F77_FUNC(dtbmv)("U", "N", "N", &n, &k, px, &incx, py, &incy);
}


void
la_vector_div (LAVector *dst_vector,
LAVector *vector)
{
    g_assert(dst_vector);
    g_assert(vector);
    g_assert(la_vector_dim(dst_vector) == la_vector_dim(vector));

    la_size n     = la_vector_dim(dst_vector);
    la_size k     = 0;
    double *px    = la_vector_ptr(vector, 0);
    la_size incx  = 1;
    double *py    = la_vector_ptr(dst_vector, 0);
    la_size incy  = 1;

    F77_FUNC(dtbsv)("U", "N", "N", &n, &k, px, &incx, py, &incy);
}


void
la_vector_acc (LAVector *dst_vector,
double    scale,
LAVector *vector)
{
    g_assert(dst_vector);
    g_assert(vector);
    g_assert(la_vector_dim(dst_vector) == la_vector_dim(vector));

    la_size n     = la_vector_dim(dst_vector);
    double  alpha = scale;
    double *px    = la_vector_ptr(vector, 0);
    la_size incx  = 1;
    double *py    = la_vector_ptr(dst_vector, 0);
    la_size incy  = 1;

    F77_FUNC(daxpy)(&n, &alpha, px, &incx, py, &incy);
}


double
la_vector_dot (LAVector *vector1,
LAVector *vector2)
{
    g_assert(vector1);
    g_assert(vector2);
    g_assert(la_vector_dim(vector1) == la_vector_dim(vector2));

    la_size n    = la_vector_dim(vector1);
    double *px   = la_vector_ptr(vector1, 0);
    la_size incx = 1;
    double *py   = la_vector_ptr(vector2, 0);
    la_size incy = 1;

    double dot = F77_FUNC(ddot)(&n, px, &incx, py, &incy);
    return dot;
}


double
la_vector_norm (LAVector *vector)
{
    g_assert(vector);

    la_size n    = la_vector_dim(vector);
    void   *px   = la_vector_ptr(vector, 0);
    la_size incx = 1;

    double norm = F77_FUNC(dnrm2)(&n, px, &incx);
    return norm;
}


double
la_vector_sum_abs (LAVector *vector)
{
    g_assert(vector);

    la_size n    = la_vector_dim(vector);
    void   *px   = la_vector_ptr(vector, 0);
    la_size incx = 1;

    double sum_abs = F77_FUNC(dasum)(&n, px, &incx);
    return sum_abs;
}


double
la_vector_max_abs (LAVector *vector)
{
    g_assert(vector);

    double max_abs = 0;
    double e;
    la_index i;

    if (la_vector_dim(vector) > 0) {
        i = la_vector_max_abs_index(vector);
        e = la_vector_get_real(vector, i);
        max_abs = fabs(e);
    }

    return max_abs;
}


la_index
la_vector_max_abs_index (LAVector *vector)
{
    g_assert(vector);
    g_return_val_if_fail(la_vector_dim(vector) > 0, -1);

    la_size n    = la_vector_dim(vector);
    void   *px   = la_vector_ptr(vector, 0);
    la_size incx = 1;

    la_index index1 = F77_FUNC(idamax)(&n, px, &incx);
    la_index index  = index - 1;
    return index;
}



