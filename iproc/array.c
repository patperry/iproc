#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include <iproc/memory.h>
#include <iproc/array.h>

static void
iproc_array_grow (iproc_array *array)
{
    assert(array);
    assert(array->n_max > 0);

    int64_t n_max = 2 * (array->n_max);
    size_t elem_size = iproc_array_elem_size(array);
    array->data = iproc_realloc(array->data, n_max * elem_size);
    array->n_max = n_max;
}

static void
iproc_array_reserve (iproc_array *array, int64_t n)
{
    assert(array);

    while (array->n_max < n) {
        iproc_array_grow(array);
    }
}

iproc_array *
iproc_array_new (size_t elem_size)
{
    assert(elem_size > 0);

    iproc_array *array = iproc_malloc(sizeof(*array));

    if (!array) return NULL;
    array->elem_size = elem_size;
    array->n = 0;
    array->n_max = 1;
    array->data = iproc_malloc(array->n_max * elem_size);

    return array;
}

void
iproc_array_free (iproc_array *array)
{
    if (array) {
        iproc_free(array->data);
        iproc_free(array);
    }
}

size_t
iproc_array_elem_size (iproc_array *array)
{
    assert(array);
    return array->elem_size;
}

void
iproc_array_set_size (iproc_array *array,
                      int64_t      n)
{
    assert(array);
    assert(n >= 0);

    int64_t nold = iproc_array_size(array);
    size_t elem_size = iproc_array_elem_size(array);

    if (n > nold) {
        iproc_array_reserve(array, n);
        memset(&(iproc_array_index(array, char, nold * elem_size)), 0,
               (n - nold) * elem_size);
    }

    array->n = n;
}

int64_t
iproc_array_size (iproc_array *array)
{
    assert(array);
    return array->n;
}

void
iproc_array_append (iproc_array *array,
                    void        *pe)
{
    assert(array);
    assert(pe);
    int64_t n = iproc_array_size(array);
    iproc_array_insert(array, n, pe);
}

void
iproc_array_prepend (iproc_array *array,
                     void        *pe)
{
    assert(array);
    assert(pe);
    iproc_array_insert(array, 0, pe);
}

void
iproc_array_insert (iproc_array *array,
                    int64_t      i,
                    void        *pe)
{
    assert(array);
    assert(pe);
    assert(i >= 0);
    int64_t n = iproc_array_size(array);
    size_t elem_size = iproc_array_elem_size(array);
    void *dst = &(iproc_array_index(array, char, i * elem_size));

    if (i < n) {
        iproc_array_set_size(array, n + 1);
        memmove(dst + elem_size, dst, (n - i) * elem_size);
    } else {
        iproc_array_set_size(array, i + 1);
    }

    memcpy(dst, pe, elem_size);
}

void
iproc_array_remove (iproc_array *array,
                    int64_t      i)
{
    assert(array);
    assert(0 <= i);
    assert(i < iproc_array_size(array));

    int64_t n = iproc_array_size(array);
    size_t elem_size = iproc_array_elem_size(array);
    void *ptr = &(iproc_array_index(array, char, i * elem_size));

    memmove(ptr, ptr + elem_size, (n - i) * elem_size);

    iproc_array_set_size(array, n - 1);
}
