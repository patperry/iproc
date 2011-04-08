#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "darray.h"


#define INITIAL_CAPACITY 1


struct darray * _darray_init (struct darray *a, size_t elt_size)
{
    assert(a);
    assert(elt_size > 0);
    
    if (_array_init(&a->array, INITIAL_CAPACITY, elt_size)) {
        a->size = 0;
        return a;
    }
    
    return NULL;
}


struct darray * darray_init_copy (struct darray *a,
                                  const struct darray *src)
{
    assert(a);
    assert(src);
    
    if (_darray_init(a, darray_elt_size(src))) {
        darray_copy(src, a);
        return a;
    }
    
    return NULL;
}


void darray_deinit (struct darray *a)
{
    assert(a);
    array_deinit(&a->array);
}


void * darray_assign (struct darray *a, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
    assert(val);
    
    darray_clear(a);
    darray_reserve(a, n);
    return darray_insert_many(a, 0, n, val);
}


void * darray_assign_array (struct darray *a, const void *ptr, ssize_t n)
{
    assert(a);
    assert(ptr || n == 0);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
    
    darray_clear(a);
    darray_reserve(a, n);
    return darray_insert_array(a, 0, ptr, n);
}


void darray_copy (const struct darray *a, struct darray *dst)
{
    assert(a);
    assert(dst);

    const void *ptr = darray_begin(a);
    ssize_t n = darray_size(a);

    darray_assign_array(dst, ptr, n);
}


void darray_swap (struct darray *a, struct darray *b)
{
    assert(a);
    assert(b);

    struct darray tmp = *a;
    *a = *b;
    *b = tmp;
}


static void darray_insert_space (struct darray *a, ssize_t i, ssize_t n)
{
    ssize_t size0 = a->size;
    ssize_t size = size0 + n;
    size_t elt_size = darray_elt_size(a);
    size_t tail_size = (size0 - i) * elt_size;

    void *src, *dst;
    
    darray_resize(a, size);
    src = darray_ptr(a, i); // compute dst after resize in case of realloc
    dst = darray_ptr(a, i + n);
    memmove(dst, src, tail_size);
}


void * darray_insert (struct darray *a, ssize_t i, const void *val)
{
    assert(a);
    assert(darray_size(a) <= darray_max_size(a) - 1);
    assert(i >= 0);
    assert(i <= darray_size(a));
    assert(val);

    return darray_insert_array(a, i, val, 1);
}


void * darray_insert_many (struct darray *a, ssize_t i, ssize_t n,
                           const void *val)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a));
    assert(val || n == 0);
    assert(n >= 0);
    assert(n <= darray_max_size(a) - darray_size(a));
    
    size_t elt_size = darray_elt_size(a);
    void *dst, *end;
    
    darray_insert_space(a, i, n);
    dst = darray_ptr(a, i);
    end = dst + i * elt_size;

    for (; dst < end; dst += elt_size) {
        memcpy(dst, val, elt_size);
    }
    
    return (void *)val + elt_size;
}


void * darray_insert_array (struct darray *a, ssize_t i, const void *ptr,
                            ssize_t n)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a));
    assert(ptr || n == 0);
    assert(n >= 0);
    assert(n <= darray_max_size(a) - darray_size(a));

    void *dst;
    size_t len = n * darray_elt_size(a);
    
    darray_insert_space(a, i, n);
    dst = darray_ptr(a, i);

    memcpy(dst, ptr, len);
    
    return (void *)ptr + len;
}

void * darray_push_back (struct darray *a, const void *val)
{
    assert(a);
    assert(darray_size(a) <= darray_max_size(a) - 1);
    assert(val);
    
    return darray_insert(a, darray_size(a), val);
}

void darray_erase (struct darray *a, ssize_t i)
{
    assert(a);
    assert(i >= 0);
    assert(i < darray_size(a));
    
    darray_erase_range(a, i, 1);
}

void darray_erase_range (struct darray *a, ssize_t i, ssize_t n)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a) - n);
    assert(n >= 0);

    ssize_t size = darray_size(a) - n;
    void *dst = darray_ptr(a, i);
    void *src = darray_ptr(a, i + n);
    void *end = darray_end(a);
    
    memmove(dst, src, end - src);
    darray_resize(a, size);
}

void darray_pop_back (struct darray *a)
{
    assert(a);
    assert(!darray_empty(a));
    
    darray_erase(a, darray_size(a) - 1);
}

void darray_clear (struct darray *a)
{
    assert(a);
    a->size = 0;
}

static void darray_grow (struct darray *a)
{
    assert(a);

    size_t nmax = darray_max_size(a);
    size_t n0 = darray_capacity(a);
    size_t inc = (n0 >> 1) + 1;
    size_t n = (n0 <= nmax - inc) ? n0 + inc : nmax;
    
    if (n != n0) {
        array_realloc(&a->array, n);
    }
}

void darray_reserve (struct darray *a, ssize_t n)
{
    assert(a);
    assert(n <= darray_max_size(a));
    
    while (darray_capacity(a) < n) {
        darray_grow(a);
    }
}

void darray_resize (struct darray *a, ssize_t n)
{
    assert(a);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
    
    size_t elt_size = darray_elt_size(a);
    char val[elt_size];
    memset(val, 0, elt_size);
    
    darray_resize_with(a, n, val);
}


void darray_resize_with (struct darray *a, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
    assert(val);

    ssize_t n0 = a->size;
    
    if (n > n0) {
        size_t elt_size = darray_elt_size(a);
        void *dst, *end;
        
        darray_reserve(a, n);
        dst = darray_ptr(a, n0);
        end = darray_ptr(a, n);
        
        for (; dst < end; dst += elt_size) {
            memcpy(dst, val, elt_size);
        }
    }

    a->size = n;
}


ssize_t darray_find_index (const struct darray *a, ssize_t i, ssize_t n,
                           const void *key, compare_fn compar)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a) - n);
    assert(n >= 0);

    return array_find_index(&a->array, i, n, key, compar);
}


ssize_t darray_find_last_index (const struct darray *a, ssize_t i, ssize_t n,
                                const void *key, compare_fn compar)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a) - n);
    assert(n >= 0);

    return array_find_last_index(&a->array, i, n, key, compar);
}


ssize_t darray_binary_search (const struct darray *a, ssize_t i, ssize_t n,
                              const void *key, compare_fn compar)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a) - n);
    assert(n >= 0);

    return array_binary_search(&a->array, i, n, key, compar);
}

