#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <iproc/blas-private.h>
#include <iproc/memory.h>
#include <iproc/svector.h>

static int
compare_int64 (void *px, void *py)
{
    int64_t x = *(int64_t *)px;
    int64_t y = *(int64_t *)py;

    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

int64_t
iproc_svector_find_nz (iproc_svector *svector, int64_t i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_dim(svector));
    int64_t ix = iproc_array_bsearch(svector->index, &i, compare_int64);
    return ix;
}

static void
iproc_svector_free (iproc_svector *svector)
{
    if (svector) {
        iproc_array_unref(svector->index);
        iproc_array_unref(svector->value);
        iproc_free(svector);
    }
}

iproc_svector *
iproc_svector_new (int64_t dim)
{
    assert(dim >= 0);
    iproc_svector *svector;

    if (!(dim >= 0)) return NULL;

    svector = iproc_malloc(sizeof(*svector));
    if (!svector) return NULL;

    svector->dim = dim;
    svector->index = iproc_array_new(sizeof(int64_t));
    svector->value = iproc_array_new(sizeof(double));
    iproc_refcount_init(&svector->refcount);

    if (!(svector->index && svector->value)) {
        iproc_svector_free(svector);
        svector = NULL;
    }

    return svector;
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
iproc_svector_new_copy (iproc_svector *svector)
{
    assert(svector);
    int64_t dim = iproc_svector_dim(svector);
    iproc_svector *copy = iproc_malloc(sizeof(*copy));
    copy->dim = dim;
    copy->index = iproc_array_new_copy(svector->index);
    copy->value = iproc_array_new_copy(svector->value);
    iproc_refcount_init(&copy->refcount);
    return copy;
}

void
iproc_svector_clear (iproc_svector *svector)
{
    assert(svector);
    iproc_array_set_size(svector->index, 0);
    iproc_array_set_size(svector->value, 0);
}


int64_t
iproc_svector_dim (iproc_svector *svector)
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
        value = iproc_svector_nz_val(svector, ix);
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
        iproc_array_insert(svector->index, ix, &i);
        iproc_array_insert(svector->value, ix, &value);
    } else {
        iproc_array_index(svector->value, double, ix) = value;
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
        iproc_array_insert(svector->index, ix, &i);
        iproc_array_insert(svector->value, ix, &value);
    } else {
        iproc_array_index(svector->value, double, ix) += value;
    }
}

void
iproc_svector_scale (iproc_svector *svector,
                     double         scale)
{
    assert(svector);

    int64_t  n     = iproc_svector_nnz(svector);
    double   alpha = scale;
    int64_t  incx  = 1;
    double  *px    = &(iproc_array_index(svector->value, double, 0));

    F77_FUNC(dscal)(F77_INTP(n), &alpha, px, F77_INTP(incx));
}

int64_t
iproc_svector_nnz (iproc_svector *svector)
{
    assert(svector);
    return iproc_array_size(svector->index);
}

int64_t
iproc_svector_nz (iproc_svector *svector,
                  int64_t        i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_nnz(svector));
    return iproc_array_index(svector->index, int64_t, i);
}

double
iproc_svector_nz_val (iproc_svector *svector,
                      int64_t        i)
{
    assert(svector);
    assert(0 <= i);
    assert(i < iproc_svector_nnz(svector));
    return iproc_array_index(svector->value, double, i);
}

iproc_vector_view
iproc_svector_view_nz (iproc_svector *svector)
{
    assert(svector);

    int64_t  nnz = iproc_svector_nnz(svector);
    double  *px  = &(iproc_array_index(svector->value, double, 0));
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
        e2 = iproc_svector_nz_val(svector, i);
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
        e = iproc_svector_nz_val(svector, i);
        ix = iproc_svector_nz(svector, i);
        iproc_vector_inc(dst_vector, ix, scale * e);
    }
}
