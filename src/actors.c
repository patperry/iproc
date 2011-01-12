#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include "memory.h"
#include "actors.h"

static int64_t
iproc_actors_unsafe_append_class (iproc_actors *actors,
                                  iproc_vector *traits)
{
    assert(actors);
    assert(traits);

    iproc_vector *x = iproc_vector_new_copy(traits);
    iproc_array_append(actors->class_traits, &x);
    return (iproc_array_size(actors->class_traits) - 1);
}

static void
iproc_actors_class_traits_free (iproc_array *class_traits)
{
    int64_t i, n;
    iproc_vector *x;

    if (class_traits) {
        n = iproc_array_size(class_traits);
        for (i = 0; i < n; i++) {
            x = iproc_array_index(class_traits, iproc_vector *, i);
            iproc_vector_unref(x);
        }
        iproc_array_unref(class_traits);
    }
}

static void
iproc_actors_free (iproc_actors *actors)
{
    if (actors) {
        iproc_actors_class_traits_free(actors->class_traits);
        iproc_array_unref(actors->classes);
        iproc_free(actors);
    }
}

iproc_actors *
iproc_actors_new (int64_t      size,
                 iproc_vector *traits0)
{
    assert(0 <= size);
    assert(traits0);
    iproc_actors *actors = iproc_malloc(sizeof(*actors));
    if (!actors) return NULL;

    actors->class_traits = iproc_array_new(sizeof(iproc_vector *));
    actors->classes = iproc_array_new(sizeof(int64_t));
    iproc_refcount_init(&actors->refcount);

    if (!(actors->class_traits && actors->classes)) {
        iproc_actors_free(actors);
        actors = NULL;
    }
    /* We rely on array_set growing clearing the tail to zeros.  No further
     * action is necessary since IPROC_ACTORS_DEFCLASS is zero.
     */
    iproc_array_set_size(actors->classes, size);
    iproc_actors_unsafe_append_class(actors, traits0);

    return actors;
}

iproc_actors *
iproc_actors_ref (iproc_actors *actors)
{
    if (actors) {
        iproc_refcount_get(&actors->refcount);
    }
    return actors;
}

static void
iproc_actors_release (iproc_refcount *refcount)
{
    iproc_actors *actors = container_of(refcount, iproc_actors, refcount);
    iproc_actors_free(actors);
}

void
iproc_actors_unref (iproc_actors *actors)
{
    if (!actors)
        return;

    iproc_refcount_put(&actors->refcount, iproc_actors_release);
}

int64_t
iproc_actors_append_class (iproc_actors *actors,
                           iproc_vector *traits)
{
    assert(actors);
    assert(traits);
    assert(iproc_vector_dim(traits) == iproc_actors_dim(actors));
    return iproc_actors_unsafe_append_class(actors, traits);
}


iproc_vector *
iproc_actors_class_traits (iproc_actors *actors,
                           int64_t       c)
{
    assert(actors);
    assert(0 <= c);
    assert(c < iproc_actors_nclass(actors));

    iproc_vector *x = iproc_array_index(actors->class_traits,
                                        iproc_vector *,
                                        c);
    return x;
}

int64_t
iproc_actors_nclass (iproc_actors *actors)
{
    assert(actors);
    return iproc_array_size(actors->class_traits);
}

int64_t
iproc_actors_size (iproc_actors *actors)
{
    assert(actors);
    return iproc_array_size(actors->classes);
}

int64_t
iproc_actors_dim (iproc_actors *actors)
{
    assert(actors);
    iproc_vector *x0 = iproc_actors_class_traits(actors, 0);
    int64_t n = iproc_vector_dim(x0);
    return n;
}

void
iproc_actors_set (iproc_actors *actors,
                  int64_t       i,
                  int64_t       c)
{
    assert(actors);
    assert(0 <= i);
    assert(i < iproc_actors_size(actors));
    assert(0 <= c);
    assert(c < iproc_actors_nclass(actors));

    iproc_array_set(actors->classes, i, &c);
}

int64_t
iproc_actors_class (iproc_actors *actors,
                    int64_t       i)
{
    assert(actors);
    assert(0 <= i);
    assert(i < iproc_actors_size(actors));

    int64_t c = IPROC_ACTORS_DEFCLASS;
    iproc_array *classes = actors->classes;

    if (i < iproc_array_size(classes)) {
        c = iproc_array_index(classes, int64_t, i);
    }

    return c;
}

iproc_vector *
iproc_actors_traits (iproc_actors *actors,
                     int64_t       i)
{
    assert(actors);
    assert(0 <= i);
    assert(i < iproc_actors_size(actors));

    int64_t c = iproc_actors_class(actors, i);
    iproc_vector *x = iproc_actors_class_traits(actors, c);
    return x;
}

void
iproc_actors_mul (double        alpha,
                  iproc_trans   trans,
                  iproc_actors *actors,
                  iproc_vector *x,
                  double        beta,
                  iproc_vector *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_vector_dim(y) == iproc_actors_dim(actors));

    int64_t n = iproc_actors_size(actors);
    int64_t i;
    iproc_vector *row;
    double dot, entry;

    if (beta == 0) {
        iproc_vector_set_all(y, 0.0);
    } else if (beta != 1) {
        iproc_vector_scale(y, beta);
    }

    if (trans == IPROC_TRANS_NOTRANS) {
        for (i = 0; i < n; i++) {
            row = iproc_actors_traits(actors, i);
            dot = iproc_vector_dot(row, x);
            iproc_vector_inc(y, i, alpha * dot);
        }
    } else {
        for (i = 0; i < n; i++) {
            row = iproc_actors_traits(actors, i);
            entry = iproc_vector_get(x, i);
            iproc_vector_acc(y, alpha * entry, row);
        }
    }
}


void
iproc_actors_matmul (double        alpha,
                     iproc_trans   trans,
                     iproc_actors *actors,
                     iproc_matrix *x,
                     double        beta,
                     iproc_matrix *y)
{
    assert(actors);
    assert(x);
    assert(y);
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(x) == iproc_actors_dim(actors));
    assert(trans != IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(y) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(x) == iproc_actors_size(actors));
    assert(trans == IPROC_TRANS_NOTRANS
           || iproc_matrix_nrow(y) == iproc_actors_dim(actors));
    assert(iproc_matrix_ncol(x) == iproc_matrix_ncol(y));

    int64_t m = iproc_matrix_ncol(x);
    int64_t j;

    if (beta == 0) {
        iproc_matrix_set_all(y, 0.0);
    } else if (beta != 1) {
        iproc_matrix_scale(y, beta);
    }

    for (j = 0; j < m; j++) {
        iproc_vector_view xcol = iproc_matrix_col(x, j);
        iproc_vector_view ycol = iproc_matrix_col(y, j);
        iproc_actors_mul(alpha, trans, actors, &xcol.vector, 1.0, &ycol.vector);
    }
}
