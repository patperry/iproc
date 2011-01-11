#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include <string.h>
#include "memory.h"
#include "array.h"

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
    iproc_refcount_init(&array->refcount);

    return array;
}

iproc_array *
iproc_array_new_copy (iproc_array *array)
{
    assert(array);

    size_t elem_size = iproc_array_elem_size(array);
    iproc_array *copy = iproc_array_new(elem_size);
    int64_t n = iproc_array_size(array);

    iproc_array_set_size(copy, n);
    void *src = &(iproc_array_index(array, char, 0));
    void *dst = &(iproc_array_index(copy, char, 0));
    memcpy(dst, src, n * elem_size);

    return copy;
}

static void
iproc_array_free (iproc_array *array)
{
    if (array) {
        iproc_free(array->data);
        iproc_free(array);
    }
}

iproc_array *
iproc_array_ref (iproc_array *array)
{
    if (array) {
        iproc_refcount_get(&array->refcount);
    }
    return array;
}

static void
iproc_array_release (iproc_refcount *refcount)
{
    iproc_array *array = container_of(refcount, iproc_array, refcount);
    iproc_array_free(array);
}

void
iproc_array_unref (iproc_array *array)
{
    if (!array)
        return;

    iproc_refcount_put(&array->refcount, iproc_array_release);
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

    assert(array->n <= array->n_max);
}

void
iproc_array_set (iproc_array *array,
                 int64_t      i,
                 void        *pe)
{
    assert(array);
    assert(i >= 0);
    assert(pe);

    int64_t nold = iproc_array_size(array);
    if (i >= nold) {
        iproc_array_set_size(array, i + 1);
    }

    size_t elem_size = iproc_array_elem_size(array);
    void *dst = &(iproc_array_index(array, char, i * elem_size));
    memcpy(dst, pe, elem_size);
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
    int64_t nnew = (i < n) ? (n + 1) : (i + 1);
    void *dst;

    iproc_array_set_size(array, nnew);

    /* compute dst after call to set_size, since memory gets realloced */
    dst = &(iproc_array_index(array, char, i * elem_size));

    if (i < n) {
        memmove(dst + elem_size, dst, (n - i) * elem_size);
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

int64_t
iproc_array_lfind (iproc_array *array,
                   void        *value,
                   int        (*compare) (void *, void *))
{
    assert(array);
    assert(value);
    assert(compare);

    char *ptr = &(iproc_array_index(array, char, 0));
    size_t elem_size = iproc_array_elem_size(array);
    int64_t n = iproc_array_size(array);
    int64_t i;

    for (i = 0; i < n; ptr += elem_size, i++) {
        if (compare(value, ptr) == 0)
            return i;
    }

    return ~i;
}

int64_t
iproc_array_bsearch (iproc_array *array,
                     void        *value,
                     int        (*compare) (void *, void *))
{
    assert(array);
    assert(value);
    assert(compare);

    size_t  elem_size = iproc_array_elem_size(array);
    void   *base = &(iproc_array_index(array, char, 0));
    int64_t begin = 0;
    int64_t end = iproc_array_size(array);
    int64_t i;
    void *ptr;
    int cmp;

    while (begin < end) {
        i = begin + ((end - begin) >> 1);
        ptr = base + i * elem_size;
        cmp = compare(value, ptr);

        if (cmp < 0) {         /* value < array[i] */
            end = i;
        } else if (cmp > 0) {  /* value > array[i] */
            begin = i + 1;
        } else {               /* value == array[i] */
            return i;
        }
    }

    return ~end;               /* begin == end, not found */
}
