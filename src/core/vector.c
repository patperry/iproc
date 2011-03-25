
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "blas-private.h"
#include "compare.h"
#include "hash.h"
#include "ieee754.h"
#include "memory.h"
#include "vector.h"


iproc_vector *
iproc_vector_new (int64_t dim)
{
    iproc_vector *vector = NULL;
    size_t        size   = dim * sizeof(double);

    if (!(dim >= 0)) return NULL;
    if (!(dim <= F77_INT_MAX)) return NULL;
    if (!(dim <= SIZE_MAX / sizeof(double))) return NULL;

    vector = iproc_malloc(sizeof(iproc_vector));
    vector->pdata = iproc_malloc(size);
    vector->dim = dim;
    iproc_refcount_init(&vector->refcount);

    return vector;
}

iproc_vector *
iproc_vector_new_copy (iproc_vector *vector)
{
    assert(vector);
    int64_t n = iproc_vector_dim(vector);
    iproc_vector *copy = iproc_vector_new(n);
    iproc_vector_copy(copy, vector);
    return copy;
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
        iproc_refcount_get(&vector->refcount);
    }

    return vector;
}

static void
iproc_vector_release (iproc_refcount *refcount)
{
    iproc_vector *vector = container_of(refcount, iproc_vector, refcount);
    iproc_vector_free(vector);
}

void
iproc_vector_unref (iproc_vector *vector)
{
    if (!vector)
        return;

    iproc_refcount_put(&vector->refcount, iproc_vector_release);
}

int64_t
iproc_vector_dim (iproc_vector *vector)
{
    if (!vector)
        return 0;
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

    iproc_vector_view view = {{ array, dim, {0} }};
    return view;
}


void
iproc_vector_copy (iproc_vector *dst_vector,
                   iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    f77int  n    = (f77int)iproc_vector_dim(dst_vector);
    double *px   = iproc_vector_ptr(vector, 0);
    f77int  incx = 1;
    double *py   = iproc_vector_ptr(dst_vector, 0);
    f77int  incy = 1;

    F77_FUNC(dcopy)(&n, px, &incx, py, &incy);
}


void
iproc_vector_swap (iproc_vector *vector1,
                   iproc_vector *vector2)
{
    assert(vector1);
    assert(vector2);
    assert(iproc_vector_dim(vector1) == iproc_vector_dim(vector2));

    f77int  n    = (f77int)iproc_vector_dim(vector1);
    double *px   = iproc_vector_ptr(vector1, 0);
    f77int  incx = 1;
    double *py   = iproc_vector_ptr(vector2, 0);
    f77int  incy = 1;

    F77_FUNC(dswap)(&n, px, &incx, py, &incy);
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

    f77int  n     = (f77int)iproc_vector_dim(vector);
    double  alpha = scale;
    void   *px    = iproc_vector_ptr(vector, 0);
    f77int  incx  = 1;

    F77_FUNC(dscal)(&n, &alpha, px, &incx);
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

    f77int  n     = (f77int)iproc_vector_dim(dst_vector);
    f77int  k     = 0;
    double *px    = iproc_vector_ptr(vector, 0);
    f77int  incx  = 1;
    double *py    = iproc_vector_ptr(dst_vector, 0);
    f77int  incy  = 1;

    F77_FUNC(dtbmv)("U", "N", "N", &n, &k, px, &incx, py, &incy);
}


void
iproc_vector_div (iproc_vector *dst_vector,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    f77int  n     = (f77int)iproc_vector_dim(dst_vector);
    f77int  k     = 0;
    double *px    = iproc_vector_ptr(vector, 0);
    f77int  incx  = 1;
    double *py    = iproc_vector_ptr(dst_vector, 0);
    f77int  incy  = 1;

    F77_FUNC(dtbsv)("U", "N", "N", &n, &k, px, &incx, py, &incy);
}


void
iproc_vector_acc (iproc_vector *dst_vector,
                  double        scale,
                  iproc_vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(iproc_vector_dim(dst_vector) == iproc_vector_dim(vector));

    f77int  n     = (f77int)iproc_vector_dim(dst_vector);
    double  alpha = scale;
    double *px    = iproc_vector_ptr(vector, 0);
    f77int  incx  = 1;
    double *py    = iproc_vector_ptr(dst_vector, 0);
    f77int  incy  = 1;

    F77_FUNC(daxpy)(&n, &alpha, px, &incx, py, &incy);
}


double
iproc_vector_dot (iproc_vector *vector1,
                  iproc_vector *vector2)
{
    assert(vector1);
    assert(vector2);
    assert(iproc_vector_dim(vector1) == iproc_vector_dim(vector2));

    f77int  n    = (f77int)iproc_vector_dim(vector1);
    double *px   = iproc_vector_ptr(vector1, 0);
    f77int  incx = 1;
    double *py   = iproc_vector_ptr(vector2, 0);
    f77int  incy = 1;

    double dot = F77_FUNC(ddot)(&n, px, &incx, py, &incy);
    return dot;
}


double
iproc_vector_norm (iproc_vector *vector)
{
    assert(vector);

    f77int  n    = (f77int)iproc_vector_dim(vector);
    void   *px   = iproc_vector_ptr(vector, 0);
    f77int  incx = 1;

    double norm = F77_FUNC(dnrm2)(&n, px, &incx);
    return norm;
}


double
iproc_vector_sum_abs (iproc_vector *vector)
{
    assert(vector);

    f77int  n    = (f77int)iproc_vector_dim(vector);
    void   *px   = iproc_vector_ptr(vector, 0);
    f77int  incx = 1;

    double sum_abs = F77_FUNC(dasum)(&n, px, &incx);
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

    f77int  n    = (f77int)iproc_vector_dim(vector);
    void   *px   = iproc_vector_ptr(vector, 0);
    f77int  incx = 1;

    f77int  index1 = F77_FUNC(idamax)(&n, px, &incx);
                                      
    int64_t index  = index1 - 1;
    return index;
}

int64_t
iproc_vector_max_index (iproc_vector *vector)
{
    assert(vector);
    assert(iproc_vector_dim(vector) > 0);

    int64_t n = iproc_vector_dim(vector);
    int64_t i, imax;
    double x, max = NAN;

    /* Find the first non-NaN entry of the vector */
    for (imax = 0; imax < n && isnan(max); imax++) {
        max = iproc_vector_get(vector, imax);
    }

    /* If all of the entries are NaN, define imax as 0. */
    if (imax == n)
        return 0;

    /* Otherwise, search for the largest entry in the tail of the vector */
    for (i = imax + 1; i < n; i++) {
        x = iproc_vector_get(vector, i);
        if (x > max) {
            max = x;
            imax = i;
        }
    }

    return imax;
}

double
iproc_vector_max (iproc_vector *vector)
{
    assert(vector);
    if (iproc_vector_dim(vector) == 0) {
        return -INFINITY;
    } else {
        int64_t i = iproc_vector_max_index(vector);
        return iproc_vector_get(vector, i);
    }
}

void
iproc_vector_exp (iproc_vector *vector)
{
    assert(vector);
    int64_t n = iproc_vector_dim(vector);
    int64_t i;
    double x;

    for (i = 0; i < n; i++) {
        x = iproc_vector_get(vector, i);
        iproc_vector_set(vector, i, exp(x));
    }
}

double 
iproc_vector_log_sum_exp (iproc_vector *vector)
{
    assert(vector);
    int64_t n = iproc_vector_dim(vector);

    if (n == 0)
        return -INFINITY;

    int64_t imax = iproc_vector_max_index(vector);
    double max = iproc_vector_get(vector, imax);
    double summ1 = 0.0;
    int64_t i;

    for (i = 0; i < n; i++) {
        if (i == imax)
            continue;

        summ1 += exp(iproc_vector_get(vector, i) - max);
    }

    return max + log1p(summ1);
}

void
iproc_vector_printf (iproc_vector *vector)
{
    printf("\nvector {");
    printf("\n  dim: %"PRId64"", iproc_vector_dim(vector));
    printf("\n   nz: {");

    int64_t i, n = iproc_vector_dim(vector);
    for (i = 0; i < n; i++) {
        if (iproc_vector_get(vector, i) == 0.0)
            continue;

        printf("\n         %"PRId64", %.8f", i,
               iproc_vector_get(vector, i));
    }
    printf("\n       }");
    printf("\n}\n");
}

size_t
iproc_vector_hash (iproc_vector *vector)
{
    if (!vector)
        return 0;

    size_t seed = 0;
    int64_t i, n = iproc_vector_dim(vector);

    for (i = 0; i < n; i++) {
        double v = iproc_vector_get(vector, i);
        size_t hash_value = iproc_hash_double(v);
        seed = iproc_hash_combine(seed, hash_value);
    }

    return seed;
}

int
iproc_vector_identical (iproc_vector *vector1,
                        iproc_vector *vector2)
{
    if (vector1 == vector2)
        return 1;

    int64_t n = iproc_vector_dim(vector1);

    if (iproc_vector_dim(vector2) != n)
        return 0;

    int64_t i;
    for (i = 0; i < n; i++) {
        double x1 = iproc_vector_get(vector1, i);
        double x2 = iproc_vector_get(vector2, i);

        if (!iproc_identical(x1, x2))
            return 0;
    }

    return 1;
}

int
iproc_vector_compare (void *x1, void *x2)
{
    iproc_vector *vector1 = x1;
    iproc_vector *vector2 = x2;
    int64_t n1 = iproc_vector_dim(vector1);
    int64_t n2 = iproc_vector_dim(vector2);

    if (n1 < n2) {
        return -1;
    } else if (n1 > n2) {
        return +1;
    }

    double *p1 = iproc_vector_ptr(vector1, 0);
    double *p2 = iproc_vector_ptr(vector2, 0);
    int64_t i, n = n1;

    for (i = 0; i < n; i++) {
        int cmp = iproc_double_compare(p1 + i, p2 + i);
        if (cmp != 0)
            return cmp;
    }

    return 0;
}

int
iproc_vector_ptr_compare (void *px1, void *px2)
{
    iproc_vector **pvector1 = px1;
    iproc_vector **pvector2 = px2;
    
    if (!pvector1)
        return pvector2 ? 0 : -1;
    if (!pvector2)
        return +1;
    
    return iproc_vector_compare(*pvector1, *pvector2);
}

