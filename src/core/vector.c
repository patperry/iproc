#include "port.h"

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

struct vector * vector_init (struct vector *v, ssize_t n)
{
    assert(v);
    assert(n >= 0);
    assert(n <= F77_INT_MAX);
    assert(n <= SSIZE_MAX / sizeof(double));

    if (array_init(&v->array, double, n)) {
        return v;
    }
    
    return NULL;
}


struct vector * vector_init_view (struct vector *v, double *ptr, ssize_t n)
{
    assert(v);
    assert(ptr || n == 0);
    assert(n >= 0);
    
    if (array_init_view(&v->array, double, ptr, n)) {
        return v;
    }

    return NULL;
}

struct vector * vector_init_slice (struct vector *v,
                                   const struct vector *parent,
                                   ssize_t i, ssize_t n)
{
    assert(v);
    assert(parent);
    assert(n >= 0);
    assert(0 <= i && i <= vector_dim(parent) - n);
    
    return vector_init_view(v, vector_ptr(parent, i), n);
}


struct vector * vector_init_copy (struct vector *v, const struct vector *src)
{
    assert(v);
    assert(src);
    
    if (vector_init(v, vector_dim(src))) {
        vector_copy(v, src);
        return v;
    }
    
    return NULL;
}


struct vector * vector_new (ssize_t n)
{
    struct vector *v = malloc(sizeof(*v));
    
    if (vector_init(v, n)) {
        return v;
    }
    
    free(v);
    return NULL;
}


struct vector * vector_new_copy (const struct vector *src)
{
    struct vector *v = malloc(sizeof(*v));
    
    if (vector_init_copy(v, src)) {
        return v;
    }
    
    free(v);
    return NULL;
}


void vector_deinit (struct vector *v)
{
    assert(v);
    array_deinit(&v->array);
}


void vector_free (struct vector *v)
{
    if (v) {
        vector_deinit(v);
        free(v);
    }
}


ssize_t
vector_dim (const struct vector *v)
{
    if (!v)
        return 0;
    return array_size(&v->array);
}


void
vector_fill (struct vector *vector,
                      double        value)
{
    assert(vector);

    ssize_t n   = vector_dim(vector);
    double *ptr = vector_ptr(vector, 0);
    double *end = vector_ptr(vector, n);

    while (ptr != end) {
        *ptr++ = value;
    }
}


void
vector_set_basis (struct vector *vector,
                        ssize_t       index)
{
    assert(vector);
    assert(vector_dim(vector) > 0);
    assert(0 <= index);
    assert(index < vector_dim(vector));

    vector_fill(vector, 0);
    vector_set(vector, index, 1);
}


double
vector_get (struct vector *vector,
                  ssize_t       index)
{
    assert(vector);
    assert(0 <= index);
    assert(index < vector_dim(vector));

    double *ptr = vector_ptr(vector, index);
    return *ptr;
}


void
vector_set (struct vector *vector,
                  ssize_t       index,
                  double        value)
{
    assert(vector);
    assert(0 <= index);
    assert(index < vector_dim(vector));

    double *ptr = vector_ptr(vector, index);
    *ptr = value;
}

void
vector_inc (struct vector *vector,
                  ssize_t       index,
                  double        value)
{
    assert(vector);
    assert(0 <= index);
    assert(index < vector_dim(vector));

    double *ptr = vector_ptr(vector, index);
    *ptr += value;
}

double *
vector_ptr (const struct vector *v, ssize_t i)
{
    assert(v);
    return &array_index(&v->array, double, i);
}


iproc_vector_view
vector_slice (struct vector *v, ssize_t i, ssize_t n)
{
    iproc_vector_view view;
    vector_init_slice(&view.vector, v, i, n);
    return view;
}


iproc_vector_view
iproc_vector_view_array (double  *ptr, ssize_t n)
{
    iproc_vector_view view;
    vector_init_view(&view.vector, ptr, n);
    return view;
}


void
vector_copy (struct vector *dst_vector,
                   const struct vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(vector_dim(dst_vector) == vector_dim(vector));

    f77int  n    = (f77int)vector_dim(dst_vector);
    double *px   = vector_ptr(vector, 0);
    f77int  incx = 1;
    double *py   = vector_ptr(dst_vector, 0);
    f77int  incy = 1;

    F77_FUNC(dcopy)(&n, px, &incx, py, &incy);
}


void
vector_swap (struct vector *vector1,
                   struct vector *vector2)
{
    assert(vector1);
    assert(vector2);
    assert(vector_dim(vector1) == vector_dim(vector2));

    f77int  n    = (f77int)vector_dim(vector1);
    double *px   = vector_ptr(vector1, 0);
    f77int  incx = 1;
    double *py   = vector_ptr(vector2, 0);
    f77int  incy = 1;

    F77_FUNC(dswap)(&n, px, &incx, py, &incy);
}


void
vector_swap_elems (struct vector *vector,
                         ssize_t       index1,
                         ssize_t       index2)
{
    assert(vector);
    assert(0 <= index1);
    assert(index1 < vector_dim(vector));
    assert(0 <= index2);
    assert(index2 < vector_dim(vector));

    double e1, e2;

    if (index1 != index2) {
        e1 = vector_get(vector, index1);
        e2 = vector_get(vector, index2);
        vector_set(vector, index2, e1);
        vector_set(vector, index1, e2);
    }
}


void
vector_reverse (struct vector *vector)
{
    assert(vector);

    ssize_t n = vector_dim(vector);
    ssize_t i;

    for (i = 0; i < n / 2; i++) {
        vector_swap_elems(vector, i, n - 1 - i);
    }
}


void
vector_scale (struct vector *vector,
                    double        scale)
{
    assert(vector);

    f77int  n     = (f77int)vector_dim(vector);
    double  alpha = scale;
    void   *px    = vector_ptr(vector, 0);
    f77int  incx  = 1;

    F77_FUNC(dscal)(&n, &alpha, px, &incx);
}


void
vector_shift (struct vector *vector,
                    double        shift)
{
    assert(vector);

    ssize_t n   = vector_dim(vector);
    double *ptr = vector_ptr(vector, 0);
    double *end = vector_ptr(vector, n);

    while (ptr != end) {
        *ptr++ += shift;
    }
}


void
vector_add (struct vector *dst_vector,
                  struct vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(vector_dim(dst_vector) == vector_dim(vector));

    vector_acc(dst_vector, 1, vector);
}


void
vector_sub (struct vector *dst_vector,
                  struct vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(vector_dim(dst_vector) == vector_dim(vector));

    vector_acc(dst_vector, -1, vector);
}


void
vector_mul (struct vector *dst_vector,
                  struct vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(vector_dim(dst_vector) == vector_dim(vector));

    f77int  n     = (f77int)vector_dim(dst_vector);
    f77int  k     = 0;
    double *px    = vector_ptr(vector, 0);
    f77int  incx  = 1;
    double *py    = vector_ptr(dst_vector, 0);
    f77int  incy  = 1;

    F77_FUNC(dtbmv)("U", "N", "N", &n, &k, px, &incx, py, &incy);
}


void
vector_div (struct vector *dst_vector,
                  struct vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(vector_dim(dst_vector) == vector_dim(vector));

    f77int  n     = (f77int)vector_dim(dst_vector);
    f77int  k     = 0;
    double *px    = vector_ptr(vector, 0);
    f77int  incx  = 1;
    double *py    = vector_ptr(dst_vector, 0);
    f77int  incy  = 1;

    F77_FUNC(dtbsv)("U", "N", "N", &n, &k, px, &incx, py, &incy);
}


void
vector_acc (struct vector *dst_vector,
                  double        scale,
                  struct vector *vector)
{
    assert(dst_vector);
    assert(vector);
    assert(vector_dim(dst_vector) == vector_dim(vector));

    f77int  n     = (f77int)vector_dim(dst_vector);
    double  alpha = scale;
    double *px    = vector_ptr(vector, 0);
    f77int  incx  = 1;
    double *py    = vector_ptr(dst_vector, 0);
    f77int  incy  = 1;

    F77_FUNC(daxpy)(&n, &alpha, px, &incx, py, &incy);
}


double
vector_dot (struct vector *vector1,
                  struct vector *vector2)
{
    assert(vector1);
    assert(vector2);
    assert(vector_dim(vector1) == vector_dim(vector2));

    f77int  n    = (f77int)vector_dim(vector1);
    double *px   = vector_ptr(vector1, 0);
    f77int  incx = 1;
    double *py   = vector_ptr(vector2, 0);
    f77int  incy = 1;

    double dot = F77_FUNC(ddot)(&n, px, &incx, py, &incy);
    return dot;
}


double
vector_norm (struct vector *vector)
{
    assert(vector);

    f77int  n    = (f77int)vector_dim(vector);
    void   *px   = vector_ptr(vector, 0);
    f77int  incx = 1;

    double norm = F77_FUNC(dnrm2)(&n, px, &incx);
    return norm;
}


double
vector_sum_abs (struct vector *vector)
{
    assert(vector);

    f77int  n    = (f77int)vector_dim(vector);
    void   *px   = vector_ptr(vector, 0);
    f77int  incx = 1;

    double sum_abs = F77_FUNC(dasum)(&n, px, &incx);
    return sum_abs;
}


double
vector_max_abs (struct vector *vector)
{
    assert(vector);

    double max_abs = 0;
    double e;
    ssize_t i;

    if (vector_dim(vector) > 0) {
        i = vector_max_abs_index(vector);
        e = vector_get(vector, i);
        max_abs = fabs(e);
    }

    return max_abs;
}


ssize_t
vector_max_abs_index (struct vector *vector)
{
    assert(vector);
    if (!(vector_dim(vector) > 0)) return -1;

    f77int  n    = (f77int)vector_dim(vector);
    void   *px   = vector_ptr(vector, 0);
    f77int  incx = 1;

    f77int  index1 = F77_FUNC(idamax)(&n, px, &incx);
                                      
    ssize_t index  = index1 - 1;
    return index;
}

ssize_t
vector_max_index (struct vector *vector)
{
    assert(vector);
    assert(vector_dim(vector) > 0);

    ssize_t n = vector_dim(vector);
    ssize_t i, imax;
    double x, max = NAN;

    /* Find the first non-NaN entry of the vector */
    for (imax = 0; imax < n && isnan(max); imax++) {
        max = vector_get(vector, imax);
    }

    /* If all of the entries are NaN, define imax as 0. */
    if (imax == n)
        return 0;

    /* Otherwise, search for the largest entry in the tail of the vector */
    for (i = imax + 1; i < n; i++) {
        x = vector_get(vector, i);
        if (x > max) {
            max = x;
            imax = i;
        }
    }

    return imax;
}

double
vector_max (struct vector *vector)
{
    assert(vector);
    if (vector_dim(vector) == 0) {
        return -INFINITY;
    } else {
        ssize_t i = vector_max_index(vector);
        return vector_get(vector, i);
    }
}

void
vector_exp (struct vector *vector)
{
    assert(vector);
    ssize_t n = vector_dim(vector);
    ssize_t i;
    double x;

    for (i = 0; i < n; i++) {
        x = vector_get(vector, i);
        vector_set(vector, i, exp(x));
    }
}

double 
vector_log_sum_exp (struct vector *vector)
{
    assert(vector);
    ssize_t n = vector_dim(vector);

    if (n == 0)
        return -INFINITY;

    ssize_t imax = vector_max_index(vector);
    double max = vector_get(vector, imax);
    double summ1 = 0.0;
    ssize_t i;

    for (i = 0; i < n; i++) {
        if (i == imax)
            continue;

        summ1 += exp(vector_get(vector, i) - max);
    }

    return max + log1p(summ1);
}

void
vector_printf (struct vector *vector)
{
    printf("\nvector {");
    printf("\n  dim: %"SSIZE_FMT"", vector_dim(vector));
    printf("\n   nz: {");

    ssize_t i, n = vector_dim(vector);
    for (i = 0; i < n; i++) {
        if (vector_get(vector, i) == 0.0)
            continue;

        printf("\n         %"SSIZE_FMT", %.8f", i,
               vector_get(vector, i));
    }
    printf("\n       }");
    printf("\n}\n");
}

size_t
vector_hash (struct vector *vector)
{
    if (!vector)
        return 0;

    size_t seed = 0;
    ssize_t i, n = vector_dim(vector);

    for (i = 0; i < n; i++) {
        double v = vector_get(vector, i);
        size_t hash_value = iproc_hash_double(v);
        seed = iproc_hash_combine(seed, hash_value);
    }

    return seed;
}

int
vector_identical (struct vector *vector1,
                        struct vector *vector2)
{
    if (vector1 == vector2)
        return 1;

    ssize_t n = vector_dim(vector1);

    if (vector_dim(vector2) != n)
        return 0;

    ssize_t i;
    for (i = 0; i < n; i++) {
        double x1 = vector_get(vector1, i);
        double x2 = vector_get(vector2, i);

        if (!iproc_identical(x1, x2))
            return 0;
    }

    return 1;
}

int
vector_compare (const void *x1, const void *x2)
{
    const struct vector *vector1 = x1;
    const struct vector *vector2 = x2;
    ssize_t n1 = vector_dim(vector1);
    ssize_t n2 = vector_dim(vector2);

    if (n1 < n2) {
        return -1;
    } else if (n1 > n2) {
        return +1;
    }

    double *p1 = vector_ptr(vector1, 0);
    double *p2 = vector_ptr(vector2, 0);
    ssize_t i, n = n1;

    for (i = 0; i < n; i++) {
        int cmp = double_compare(p1 + i, p2 + i);
        if (cmp != 0)
            return cmp;
    }

    return 0;
}

int
vector_ptr_compare (const void *px1, const void *px2)
{
    struct vector * const *pvector1 = px1;
    struct vector * const *pvector2 = px2;
    
    if (!pvector1)
        return pvector2 ? 0 : -1;
    if (!pvector2)
        return +1;
    
    return vector_compare(*pvector1, *pvector2);
}

