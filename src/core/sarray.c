#include "port.h"
#include <assert.h>
#include "sarray.h"


static ssize_t sarray_find_index (const struct sarray *a, ssize_t i)
{
    return darray_binary_search(&a->pattern, &i, ssize_compare);
}


struct sarray * _sarray_init (struct sarray *a, size_t elt_size)
{
    assert(a);
    assert(elt_size > 0);

    if (_darray_init(&a->array, elt_size)
        && darray_init(&a->pattern, ssize_t)) {
        return a;
    }
    
    return NULL;
}


struct sarray * sarray_init_copy (struct sarray *a, const struct sarray *src)
{
    assert(a);
    assert(src);

    if (darray_init_copy(&a->array, &src->array)
        && darray_init_copy(&a->pattern, &src->pattern)) {
        return a;
    }
    
    return NULL;
}


void sarray_deinit (struct sarray *a)
{
    assert(a);
    
    darray_deinit(&a->pattern);
    darray_deinit(&a->array);
}


struct sarray * sarray_assign_copy (struct sarray       *a,
                                    const struct sarray *src)
{
    assert(a);
    assert(src);
    
    if (darray_assign_copy(&a->array, &src->array)
        && darray_assign_copy(&a->pattern, &src->pattern)) {
        return a;
    }
    
    return NULL;
}


void * sarray_copy_to (const struct sarray *a, void *dst)
{
    assert(a);
    assert(dst);
    
    return darray_copy_to(&a->array, dst);
}


ssize_t * sarray_copy_pattern_to (const struct sarray *a, ssize_t *dst)
{
    assert(a);
    assert(dst);
    
    return darray_copy_to(&a->pattern, dst);
}


void * _sarray_index (struct sarray *a, ssize_t i)
{
    assert(a);
    assert(i >= 0);
    return _sarray_index_with(a, i, NULL);
}


void * _sarray_index_with (struct sarray *a, ssize_t i, const void *val0)
{
    assert(a);
    assert(i >= 0);
    
    ssize_t index = sarray_find_index(a, i);
    
    if (index < 0) {
        index = ~index;
        darray_insert(&a->array, index, val0);
        darray_insert(&a->pattern, index, &i);
    }
    
    return darray_ptr(&a->array, index);
}


void * sarray_find (const struct sarray *a, ssize_t i)
{
    assert(a);
    assert(i >= 0);

    return sarray_find_with(a, i, NULL);
}


void * sarray_find_with (const struct sarray *a, ssize_t i, const void *val0)
{
    assert(a);
    assert(i >= 0);

    ssize_t index = sarray_find_index(a, i);
    
    if (index >= 0) {
        return darray_ptr(&a->array, index);
    }
    
    return (void *)val0;
}


void sarray_replace (struct sarray *a, ssize_t i, void *src)
{
    assert(a);
    assert(i >= 0);
    
    ssize_t index = sarray_find_index(a, i);
    
    if (index >= 0) { // key already exists
        if (src) {
            darray_set(&a->array, index, src);
        } else {
            darray_erase(&a->array, index);
            darray_erase(&a->pattern, index);
        }
    } else if (src) { // key doesn't exist, src is non-NULL
        index = ~index;
        darray_insert(&a->array, index, src);
        darray_insert(&a->pattern, index, &i);
    }
}


void sarray_erase (struct sarray *a, ssize_t i)
{
    assert(a);
    assert(i >= 0);

    sarray_replace(a, i, NULL);
}


void * sarray_scatter (const struct sarray *a, void *dst)
{
    assert(a);
    assert(dst || sarray_size(a) == 0);
    
    size_t elt_size = sarray_elt_size(a);
    const void *src = sarray_begin(a);
    const void *end = sarray_end(a);
    const ssize_t *pat_ptr = sarray_pattern_begin(a);
    ssize_t i = -1;
    
    for (; src < end; src += elt_size, pat_ptr++) {
        i = *pat_ptr;
        memcpy(dst + i * elt_size, src, elt_size);
    }
    
    return dst + (i + 1) * elt_size;
}


void * sarray_gather (struct sarray *a, const void *src)
{
    assert(a);
    assert(src || sarray_size(a) == 0);
    
    size_t elt_size = sarray_elt_size(a);
    void *dst = sarray_begin(a);
    const void *end = sarray_end(a);
    const ssize_t *pat_ptr = sarray_pattern_begin(a);
    ssize_t i = -1;
    
    for (; dst < end; dst += elt_size, pat_ptr++) {
        i = *pat_ptr;
        memcpy(dst, src + i * elt_size, elt_size);
    }
    
    return (void *)src + (i + 1) * elt_size;
}

