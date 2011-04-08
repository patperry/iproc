#include "port.h"

#include <assert.h>
#include <string.h>
#include "util.h"

void * copy_to (const void *src, ssize_t size, void *dst, size_t elt_size)
{
    assert(src || size == 0);
    assert(size >= 0);
    assert(dst || size == 0);
    assert(elt_size > 0);
    
    size_t nbytes = size * elt_size;
    memmove(dst, src, nbytes); // not memcpy; be careful of aliasing
    return dst + elt_size;
}


ssize_t find_index (const void *begin, ssize_t size, const void *key,
                    compare_fn compar, size_t elt_size)
{
    assert(size >= 0);
    assert(compar);    
    assert(elt_size > 0);
    
    const void *end = begin + size * elt_size;
    const void *ptr;
    
    for (ptr = begin; ptr < end; ptr += elt_size) {
        if (compar(ptr, key) == 0)
            return (ptr - begin) / elt_size;
    }

    return -1;
}


ssize_t find_last_index (const void *begin, ssize_t size, const void *key,
                         compare_fn compar, size_t elt_size)
{
    assert(size >= 0);
    assert(compar);    
    assert(elt_size > 0);
    
    const void *end = begin + size * elt_size;
    const void *ptr;
    
    for (ptr = end; ptr > begin;) {
        ptr -= elt_size;
        if (compar(ptr, key) == 0)
            return (ptr - begin) / elt_size;
    }
    
    return -1;
}


ssize_t binary_search (const void *begin, ssize_t size, const void *key,
                       compare_fn compar, size_t elt_size)
{
    assert(size >= 0);
    assert(compar);
    assert(elt_size > 0);    
    
    ssize_t left = 0;
    ssize_t right = size;
    ssize_t i;
    const void *ptr;
    int cmp;
    
    while (left < right) {
        i = left + ((right - left) >> 1);
        ptr = begin + i * elt_size;
        cmp = compar(ptr, key);

        if (cmp < 0) {        // array[i] < key
            left = i + 1;
        } else if (cmp > 0) { // array[i] > key
            right = i;
        } else {              // array[i] = key
            return i;
        }
    }
    
    return ~right;             // left == right, not found
}
