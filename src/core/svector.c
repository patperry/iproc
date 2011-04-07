#include "port.h"

#include <assert.h>
#include <stdio.h>
#include <inttypes.h>
#include "blas-private.h"
#include "compare.h"
#include "memory.h"
#include "svector.h"


int64_t
iproc_svector_find_nz (iproc_svector *svector, int64_t i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_dim(svector));
    ssize_t ix = darray_bsearch(&svector->index, &i, int64_compare);
    return ix;
}

static void
iproc_svector_free (iproc_svector *svector)
{
    if (svector) {
        darray_deinit(&svector->index);
        darray_deinit(&svector->value);
        iproc_free(svector);
    }
}

iproc_svector *
iproc_svector_new (int64_t dim)
{
    assert(dim >= 0);

    iproc_svector *svector = iproc_calloc(1, sizeof(*svector));
    
    if (svector
        && darray_init(&svector->index, int64_t)
        && darray_init(&svector->value, double)
        && iproc_refcount_init(&svector->refcount)) {
        svector->dim = dim;
        return svector;
    }

    iproc_svector_free(svector);
    return NULL;
}

iproc_svector *
iproc_svector_ref (iproc_svector *svector)
{
    if (svector) {
        iproc_refcount_get(&svector->refcount);
    }
    return svector;
}

static void
iproc_svector_release (iproc_refcount *refcount)
{
    iproc_svector *svector = container_of(refcount, iproc_svector, refcount);
    iproc_svector_free(svector);
}

void
iproc_svector_unref (iproc_svector *svector)
{
    if (!svector)
        return;

    iproc_refcount_put(&svector->refcount, iproc_svector_release);
}

iproc_svector *
iproc_svector_new_copy (const iproc_svector *svector)
{
    assert(svector);
    int64_t dim = iproc_svector_dim(svector);
    iproc_svector *copy = iproc_calloc(1, sizeof(*copy));
    
    if (copy
        && darray_init_copy(&copy->index, &svector->index)
        && darray_init_copy(&copy->value, &svector->value)
        && iproc_refcount_init(&copy->refcount)) {
        copy->dim = dim;
        return copy;
    }

    iproc_free(copy);
    return NULL;
}

void
iproc_svector_clear (iproc_svector *svector)
{
    assert(svector);
    darray_clear(&svector->index);
    darray_clear(&svector->value);
}


int64_t
iproc_svector_dim (const iproc_svector *svector)
{
    assert(svector);
    return svector->dim;
}

double
iproc_svector_get (iproc_svector *svector,
                   int64_t        i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_dim(svector));

    double value = 0.0;
    int64_t ix = iproc_svector_find_nz(svector, i);

    if (ix >= 0) {
        value = iproc_svector_nz_get(svector, ix);
    }

    return value;
}

void
iproc_svector_set (iproc_svector *svector,
                   int64_t        i,
                   double         value)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_dim(svector));

    int64_t ix = iproc_svector_find_nz(svector, i);

    if (ix < 0) {
        ix = ~ix;
        darray_insert(&svector->index, ix, &i);
        darray_insert(&svector->value, ix, &value);
    } else {
        darray_index(&svector->value, double, ix) = value;
    }
}

void
iproc_svector_inc (iproc_svector *svector,
                   int64_t        i,
                   double         value)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_dim(svector));

    int64_t ix = iproc_svector_find_nz(svector, i);

    if (ix < 0) {
        ix = ~ix;
        darray_insert(&svector->index, ix, &i);
        darray_insert(&svector->value, ix, &value);
    } else {
        darray_index(&svector->value, double, ix) += value;
    }
}

void
iproc_svector_scale (iproc_svector *svector,
                     double         scale)
{
    assert(svector);

    f77int   n     = (f77int)iproc_svector_nnz(svector);
    double   alpha = scale;
    f77int   incx  = 1;
    double  *px    = darray_begin(&svector->value);

    F77_FUNC(dscal)(&n, &alpha, px, &incx);
}

int64_t
iproc_svector_nnz (iproc_svector *svector)
{
    assert(svector);
    return darray_size(&svector->index);
}

int64_t
iproc_svector_nz (iproc_svector *svector,
                  int64_t        i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_nnz(svector));
    return darray_index(&svector->index, int64_t, i);
}

double
iproc_svector_nz_get (iproc_svector *svector,
                      int64_t        i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_nnz(svector));
    return darray_index(&svector->value, double, i);
}

void
iproc_svector_nz_set (iproc_svector *svector,
                      int64_t        i,
                      double         value)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_nnz(svector));
    darray_index(&svector->value, double, i) = value;
}

void
iproc_svector_nz_inc (iproc_svector *svector,
                      int64_t        i,
                      double         inc)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_nnz(svector));
    darray_index(&svector->value, double, i) += inc;
}

iproc_vector_view
iproc_svector_view_nz (iproc_svector *svector)
{
    assert(svector);

    int64_t  nnz = iproc_svector_nnz(svector);
    double  *px  = &darray_index(&svector->value, double, 0);
    return iproc_vector_view_array(px, nnz);
}

double
iproc_vector_sdot (iproc_vector  *vector,
                   iproc_svector *svector)
{
    assert(vector);
    assert(svector);
    assert(iproc_vector_dim(vector) == iproc_svector_dim(svector));

    int64_t i, ix, n = iproc_svector_nnz(svector);
    double e1, e2, dot = 0.0;

    for (i = 0; i < n; i++) {
        ix = iproc_svector_nz(svector, i);
        e1 = iproc_vector_get(vector, ix);
        e2 = iproc_svector_nz_get(svector, i);
        dot += e1 * e2;
    }
    return dot;
}

void
iproc_vector_sacc (iproc_vector  *dst_vector,
                   double         scale,
                   iproc_svector *svector)
{
    assert(dst_vector);
    assert(svector);
    assert(iproc_vector_dim(dst_vector) == iproc_svector_dim(svector));

    int64_t i, ix, n = iproc_svector_nnz(svector);
    double e;

    for (i = 0; i < n; i++) {
        e = iproc_svector_nz_get(svector, i);
        ix = iproc_svector_nz(svector, i);
        iproc_vector_inc(dst_vector, ix, scale * e);
    }
}

double
iproc_svector_sdot (iproc_svector *svector1,
                    iproc_svector *svector2)
{
    assert(svector1);
    assert(svector2);
    assert(iproc_svector_dim(svector1) == iproc_svector_dim(svector2));

    int64_t n1 = iproc_svector_nnz(svector1);
    int64_t n2 = iproc_svector_nnz(svector2);

    if (n1 == 0 || n2 == 0)
        return 0.0;

    int64_t i1 = 0;
    int64_t i2 = 0;
    int64_t ix1, ix2;
    double dot = 0.0, val1, val2;

    while (i1 < n1 && i2 < n2) {
        ix1 = iproc_svector_nz(svector1, i1);
        ix2 = iproc_svector_nz(svector2, i2);

        if (ix1 == ix2) {
            val1 = iproc_svector_nz_get(svector1, i1);
            val2 = iproc_svector_nz_get(svector2, i2);
            dot += val1 * val2;
            i1++;
            i2++;
        } else if (ix1 < ix2) {
            i1++;
        } else { /* ix1 > ix2 */
            i2++;
        }
    }
    
    return dot;
}


void
iproc_svector_sacc (iproc_svector *dst_svector,
                    double         scale,
                    iproc_svector *svector)
{
    assert(dst_svector);
    assert(svector);
    assert(iproc_svector_dim(dst_svector) == iproc_svector_dim(svector));

    int64_t inz, nnz = iproc_svector_nnz(svector);
    for (inz = 0; inz < nnz; inz++) {
        int64_t i = iproc_svector_nz(svector, inz);
        double val = iproc_svector_nz_get(svector, inz);
        iproc_svector_inc(dst_svector, i, scale * val);
    }
}


void
iproc_svector_printf (iproc_svector *svector)
{
    printf("\nsvector {");
    printf("\n  dim: %"PRId64"", iproc_svector_dim(svector));
    printf("\n   nz: {");

    int64_t i, n = iproc_svector_nnz(svector);
    for (i = 0; i < n; i++) {
        printf("\n         %"PRId64", %.8f",
               iproc_svector_nz(svector, i),
               iproc_svector_nz_get(svector, i));
    }
    printf("\n       }");
    printf("\n}\n");
}

void
iproc_svector_copy (iproc_svector *dst_svector,
                    iproc_svector *src_svector)
{
    assert(dst_svector);
    assert(src_svector);
    assert(iproc_svector_dim(dst_svector) == iproc_svector_dim(src_svector));

    int64_t nnz = iproc_svector_nnz(src_svector);

    darray_resize(&dst_svector->index, nnz);
    int64_t *idst = &darray_index(&dst_svector->index, int64_t, 0);
    int64_t *isrc = &darray_index(&src_svector->index, int64_t, 0);
    memcpy(idst, isrc, nnz * sizeof(int64_t));

    darray_resize(&dst_svector->value, nnz);
    double *vdst = &darray_index(&dst_svector->value, double, 0);
    double *vsrc = &darray_index(&src_svector->value, double, 0);
    memcpy(vdst, vsrc, nnz * sizeof(double));
}
