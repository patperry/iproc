#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "array.h"


struct array * _array_init (struct array *a, ssize_t size, size_t elt_size)
{
    assert(a);
    assert(size >= 0);
    assert(elt_size > 0);
    assert(size <= SSIZE_MAX / elt_size);
    
    a->data = malloc(size * elt_size);
    a->size = size;
    a->elt_size = elt_size;
    a->owner = true;
    
    if (!a->data)
        a = NULL;
    
    return a;
}


struct array * _array_init_view (struct array *a, const void *data,
                                 ssize_t size, size_t elt_size)
{
    assert(a);
    assert(data || size == 0);
    assert(elt_size > 0);
    assert(size <= SSIZE_MAX / elt_size);

    a->data = (void *)data;
    a->size = size;
    a->elt_size = elt_size;
    a->owner = false;
    
    return a;
}


struct array * array_init_slice (struct array *a, const struct array *parent,
                                 ssize_t i, ssize_t n)
{
    assert(a);
    assert(parent);
    assert(n >= 0);
    assert(0 <= i && i <= array_size(parent) - n);
    
    return _array_init_view(a, array_ptr(parent, i), n, parent->elt_size);
}



struct array * array_init_copy (struct array *a, const struct array *src)
{
    assert(a);
    assert(src);
    
    if (_array_init(a, array_size(src), array_elt_size(src))) {
        return array_assign_copy(a, src);
    }
    
    return NULL;
}


void array_deinit (struct array *a)
{
    assert(a);
    
    if (array_owner(a)) {
        free(a->data);
    }
}


void array_swap (struct array *a, ssize_t i, ssize_t j)
{
    assert(a);
    assert(0 <= i && i < array_size(a));
    assert(0 <= j && j < array_size(a));
    assert(i != j);
    
    memory_swap(array_ptr(a, i), array_ptr(a, j), array_elt_size(a));
}


void array_reverse (struct array *a)
{
    assert(a);
    memory_reverse(array_begin(a), array_size(a), array_elt_size(a));
}


struct array * array_resize (struct array *a, ssize_t n)
{
    assert(a);
    assert(array_owner(a));
    assert(n >= 0);
    
    size_t nbytes = n * a->elt_size;
    void *data = realloc(a->data, nbytes);
    
    if (data) {
        a->data = data;
        a->size = n;
        return a;
    }
    
    return NULL;
}



struct array * array_assign_array (struct array *a, const void *ptr)
{
    assert(a);
    assert(ptr || array_size(a) == 0);

    memory_copy_to(ptr, array_size(a), array_begin(a), array_elt_size(a));
    return a;
}


struct array * array_assign_copy (struct array *a, const struct array *src)
{
    assert(a);
    assert(src);
    assert(array_elt_size(a) == array_elt_size(src));
    assert(array_size(a) == array_size(src));
    
    array_copy_to(src, array_begin(a));
    return a;
}


void * array_copy_to (const struct array *a, void *dst)
{
    assert(a);
    assert(dst || array_size(a) == 0);
    
    return memory_copy_to(array_begin(a), array_size(a), dst, array_elt_size(a));
}


void array_fill (struct array *a, const void *val)
{
    assert(a);
    array_fill_range(a, 0, array_size(a), val);
}


void array_fill_range (struct array *a, ssize_t i, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(0 <= i && i <= array_size(a) - n);

    memory_fill(array_ptr(a, i), n, val, array_elt_size(a));
}


ssize_t array_search (const struct array *a, const void *key,
                          compare_fn compar)
{
    assert(a);
    assert(compar);
    
    return forward_search(array_begin(a), array_size(a), key, compar,
                      array_elt_size(a));
}


ssize_t array_reverse_search (const struct array *a, const void *key,
                               compare_fn compar)
{
    assert(a);
    assert(compar);
    
    return reverse_search(array_begin(a), array_size(a), key, compar,
                           array_elt_size(a));
}


ssize_t array_binary_search (const struct array *a, const void *key,
                             compare_fn compar)
                             
{
    assert(a);
    assert(compar);

    return binary_search(array_begin(a), array_size(a), key, compar,
                         array_elt_size(a));
}
