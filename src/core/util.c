#include "port.h"

#include <assert.h>
#include <string.h>
#include "util.h"

void * memory_fill (void *begin, ssize_t size, const void *val, size_t elt_size)
{
    assert(begin || size == 0);
    assert(elt_size > 0);
    
    size_t nbytes = size * elt_size;
    void *end = begin + nbytes;
    void *dst;

    if (!val) {
        memset(begin, 0, nbytes);
    } else {
        assert(!(begin <= val && val < end));

        for (dst = begin; dst < end; dst += elt_size) {
            memcpy(dst, val, elt_size);
        }
    }
    
    return end;
}


void * memory_copy_to (const void *src, ssize_t size, void *dst, size_t elt_size)
{
    assert(src || size == 0);
    assert(size >= 0);
    assert(dst || size == 0);
    assert(elt_size > 0);
    
    size_t nbytes = size * elt_size;
    memmove(dst, src, nbytes); // not memcpy; be careful of aliasing
    return dst + elt_size;
}


void memory_swap (void *val1, void *val2, size_t elt_size)
{
    assert(val1);
    assert(val2);
    assert(elt_size > 0);
    assert(val1 != val2);
    
    char tmp[elt_size];

    memcpy(tmp,  val1, elt_size);
    memcpy(val1, val2, elt_size);
    memcpy(val2, tmp,  elt_size);
}


void * memory_reverse (void *begin, ssize_t size, size_t elt_size)
{
    assert(begin || size == 0);
    assert(size >= 0);
    assert(elt_size > 0);
    
    ssize_t n = size;
    ssize_t i;
    
    for (i = 0; i < n / 2; i++) {
        memory_swap(begin + i * elt_size, begin + (n - i - 1) * elt_size,
                    elt_size);
    }
    
    return begin + n * elt_size;
}


ssize_t forward_search (const void *begin, ssize_t size, const void *key,
                        equal_fn equal, size_t elt_size)
{
    assert(size >= 0);
    assert(equal);    
    assert(elt_size > 0);
    
    const void *end = begin + size * elt_size;
    const void *ptr;
    
    for (ptr = begin; ptr < end; ptr += elt_size) {
        if (equal(ptr, key))
            return (ptr - begin) / elt_size;
    }

    return -1;
}


ssize_t reverse_search (const void *begin, ssize_t size, const void *key,
                        equal_fn equal, size_t elt_size)
{
    assert(size >= 0);
    assert(equal);    
    assert(elt_size > 0);
    
    const void *end = begin + size * elt_size;
    const void *ptr;
    
    for (ptr = end; ptr > begin;) {
        ptr -= elt_size;
        if (equal(ptr, key))
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
