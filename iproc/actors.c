#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <iproc/memory.h>
#include <iproc/actors.h>

static int64_t
iproc_actors_unsafe_append_class (iproc_actors *actors,
                                  iproc_vector *vector)
{
    assert(actors);
    assert(vector);

    int64_t n = iproc_vector_dim(vector);
    iproc_vector *x = iproc_vector_new(n);
    iproc_vector_copy(x, vector);
    iproc_array_append(actors->vector, &x);
    return (iproc_array_size(actors->vector) - 1);
}

iproc_actors *
iproc_actors_new (int64_t      size,
                 iproc_vector *defvector)
{
    assert(0 <= size);
    assert(defvector);
    iproc_actors *actors = iproc_malloc(sizeof(*actors));
    if (!actors) return NULL;

    actors->vector = iproc_array_new(sizeof(iproc_vector *));
    actors->class = iproc_array_new(sizeof(int64_t));

    if (!(actors->vector && actors->class)) {
        iproc_actors_free(actors);
        actors = NULL;
    }
    /* We rely on array_set growing clearing the tail to zeros.  No further
     * action is necessary since IPROC_ACTORS_DEFCLASS is zero.
     */
    iproc_array_set_size(actors->class, size);
    iproc_actors_unsafe_append_class(actors, defvector);

    return actors;
}

static void
iproc_actors_vector_free (iproc_array *vector)
{
    int64_t i, n;
    iproc_vector *x;

    if (vector) {
        n = iproc_array_size(vector);
        for (i = 0; i < n; i++) {
            x = iproc_array_index(vector, iproc_vector *, i);
            iproc_vector_free(x);
        }
        iproc_array_unref(vector);
    }
}

void
iproc_actors_free (iproc_actors *actors)
{
    if (actors) {
        iproc_actors_vector_free(actors->vector);
        iproc_array_unref(actors->class);
        iproc_free(actors);
    }
}

int64_t
iproc_actors_append_class (iproc_actors *actors,
                           iproc_vector *vector)
{
    assert(actors);
    assert(vector);
    assert(iproc_vector_dim(vector) == iproc_actors_dim(actors));
    return iproc_actors_unsafe_append_class(actors, vector);
}


iproc_vector *
iproc_actors_class_vector (iproc_actors *actors,
                           int64_t       c)
{
    assert(actors);
    assert(0 <= c);
    assert(c < iproc_actors_nclass(actors));

    iproc_vector *x = iproc_array_index(actors->vector,
                                        iproc_vector *,
                                        0);
    return x;
}

int64_t
iproc_actors_nclass (iproc_actors *actors)
{
    assert(actors);
    return iproc_array_size(actors->vector);
}

int64_t
iproc_actors_size (iproc_actors *actors)
{
    assert(actors);
    return iproc_array_size(actors->class);
}

int64_t
iproc_actors_dim (iproc_actors *actors)
{
    assert(actors);
    iproc_vector *x0 = iproc_actors_class_vector(actors, 0);
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

    iproc_array_set(actors->class, i, &c);
}

int64_t
iproc_actors_class (iproc_actors *actors,
                    int64_t       i)
{
    assert(actors);
    assert(0 <= i);
    assert(i < iproc_actors_size(actors));

    int64_t c = IPROC_ACTORS_DEFCLASS;
    iproc_array *class = actors->class;

    if (i < iproc_array_size(class)) {
        c = iproc_array_index(class, int64_t, i);
    }

    return c;
}

iproc_vector *
iproc_actors_vector (iproc_actors *actors,
                     int64_t       i)
{
    assert(actors);
    assert(0 <= i);
    assert(i < iproc_actors_size(actors));

    int64_t c = iproc_actors_class(actors, i);
    iproc_vector *x = iproc_actors_class_vector(actors, c);
    return x;
}
