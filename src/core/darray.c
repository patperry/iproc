#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "util.h"
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
        return darray_assign_copy(a, src);
    }
    
    return NULL;
}


void darray_deinit (struct darray *a)
{
    assert(a);
    array_deinit(&a->array);
}


struct darray * darray_assign_repeat (struct darray *a, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
    assert(val);
    
    if (darray_reserve(a, n)) {
        darray_clear(a);        
        darray_insert_repeat(a, 0, n, val);
        return a;
    }

    return NULL;
}


struct darray * darray_assign (struct darray *a, const void *ptr, ssize_t n)
{
    assert(a);
    assert(ptr || n == 0);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
    
    if (darray_reserve(a, n)) {
        darray_clear(a);
        darray_insert_all(a, 0, ptr, n);
        return a;
    }

    return NULL;
}


struct darray * darray_assign_copy (struct darray *a, const struct darray *src)
{
    assert(a);
    assert(src);
    
    return darray_assign(a, darray_begin(src), darray_size(src));
}


void * darray_copy_to (const struct darray *a, void *dst)
{
    assert(a);
    assert(dst || darray_size(a) == 0);
    
    return memory_copy_to(darray_begin(a), darray_size(a), dst, darray_elt_size(a));
}


void darray_fill (struct darray *a, const void *val)
{
    assert(a);
    darray_fill_range(a, 0, darray_size(a), val);
}


void darray_fill_range (struct darray *a, ssize_t i, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(0 <= i && i <= darray_size(a) - n);
    
    memory_fill(darray_ptr(a, i), n, val, darray_elt_size(a));
}


static bool darray_insert_space (struct darray *a, ssize_t i, ssize_t n)
{
    ssize_t size0 = a->size;
    ssize_t size = size0 + n;
    size_t elt_size = darray_elt_size(a);
    size_t tail_size = (size0 - i) * elt_size;

    void *src, *dst;
    
    if (darray_resize(a, size)) {
        src = darray_ptr(a, i); // compute dst after resize in case of realloc
        dst = darray_ptr(a, i + n);
        memmove(dst, src, tail_size);
        return true;
    }

    return false;
}


void * darray_insert (struct darray *a, ssize_t i, const void *val)
{
    assert(a);
    assert(darray_size(a) <= darray_max_size(a) - 1);
    assert(i >= 0);
    assert(i <= darray_size(a));

    return darray_insert_repeat(a, i, 1, val);
}


void * darray_insert_repeat (struct darray *a, ssize_t i, ssize_t n,
                             const void *val)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a));
    assert(n >= 0);
    assert(n <= darray_max_size(a) - darray_size(a));
    
    if (darray_insert_space(a, i, n)) {
        array_fill_range(&a->array, i, n, val);
        return (void *)val + darray_elt_size(a);
    }

    return NULL;
}


void * darray_insert_all (struct darray *a, ssize_t i, const void *ptr,
                          ssize_t n)
{
    assert(a);
    assert(i >= 0);
    assert(i <= darray_size(a));
    assert(ptr || n == 0);
    assert(n >= 0);
    assert(n <= darray_max_size(a) - darray_size(a));

    if (darray_insert_space(a, i, n)) {
        return array_set_range(&a->array, i, n, ptr);
    }
    return NULL;
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


static struct darray * darray_grow (struct darray *a)
{
    assert(a);

    size_t nmax = darray_max_size(a);
    size_t n0 = darray_capacity(a);
    size_t inc = (n0 >> 1) + 1;
    size_t n = (n0 <= nmax - inc) ? n0 + inc : nmax;
    
    if (n != n0 && _array_reinit(&a->array, n, darray_elt_size(a))) {
        return a;
    }
    
    return NULL;
}


struct darray * darray_reserve (struct darray *a, ssize_t n)
{
    assert(a);
    assert(n <= darray_max_size(a));
    
    while (darray_capacity(a) < n) {
        if (!darray_grow(a))
            return NULL;
    }
    
    return a;
}

struct darray * darray_resize (struct darray *a, ssize_t n)
{
    assert(a);
    assert(n >= 0);
    assert(n <= darray_max_size(a));
        
    return darray_resize_with(a, n, NULL);
}


struct darray * darray_resize_with (struct darray *a, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(n <= darray_max_size(a));

    ssize_t n0 = a->size;
    
    if (n > n0) {
        if (!darray_reserve(a, n))
            return NULL;

        array_fill_range(&a->array, n0, n - n0, val);
    }

    a->size = n;
    return a;
}


ssize_t darray_search (const struct darray *a, const void *key,
                           compare_fn compar)
{
    assert(a);
    assert(compar);
    
    return forward_search(darray_begin(a), darray_size(a), key, compar,
                          darray_elt_size(a));
}


ssize_t darray_reverse_search (const struct darray *a, const void *key,
                                compare_fn compar)
{
    assert(a);
    assert(compar);
    
    return reverse_search(darray_begin(a), darray_size(a), key, compar,
                          darray_elt_size(a));
}


ssize_t darray_binary_search (const struct darray *a, const void *key,
                              compare_fn compar)

{
    assert(a);
    assert(compar);
    
    return binary_search(darray_begin(a), darray_size(a), key, compar,
                         darray_elt_size(a));
}

void darray_swap (struct darray *a, ssize_t i, ssize_t j)
{
    assert(a);
    assert(0 <= i && i < darray_size(a));
    assert(0 <= j && j < darray_size(a));
    assert(i != j);
    
    memory_swap(darray_ptr(a, i), darray_ptr(a, j), darray_elt_size(a));
}


void darray_reverse (struct darray *a)
{
    assert(a);
    memory_reverse(darray_begin(a), darray_size(a), darray_elt_size(a));
}

