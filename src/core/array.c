#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "array.h"


struct array * _array_init (struct array *a, ssize_t size, ssize_t elt_size)
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
                                 ssize_t size, ssize_t elt_size)
{
    assert(a);
    assert(data || size == 0);
    assert(size >= 0);
    assert(elt_size > 0);

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
    assert(i >= 0);
    assert(n >= 0);
    assert(i <= array_size(parent) - n);
    
    return _array_init_view(a, array_ptr(parent, i), n, parent->elt_size);
}



struct array * array_init_copy (struct array *a, const struct array *src)
{
    assert(a);
    assert(src);
    
    if (_array_init(a, array_size(src), _array_elt_size(src))) {
        array_copy(src, a);
        return a;
    }
    
    return NULL;
}


void array_deinit (struct array *a)
{
    assert(a);
    
    if (array_owner(a))
        free(a->data);
}


bool array_realloc (struct array *a, ssize_t n)
{
    assert(a);
    assert(array_owner(a));
    assert(n >= 0);
    
    size_t nbytes = n * a->elt_size;
    void *data = realloc(a->data, nbytes);
    
    if (data) {
        a->data = data;
        a->size = n;
    }
    
    return (data != NULL);
}



void * array_assign (struct array *a, ssize_t n, const void *val)
{
    assert(a);
    assert(n >= 0);
    assert(val || n == 0);

    char val0[a->elt_size];
    
    if (n < a->size)
        memset(val0, 0, a->elt_size);

    return array_assign_with(a, n, val, val0);
}


void * array_assign_with (struct array *a, ssize_t n, const void *val,
                          const void *val0)
{
    assert(a);
    assert(n >= 0);
    assert(val || n == 0);
    assert(val0 || n >= array_size(a));
    
    ssize_t size = a->size;
    ssize_t nsize = MIN(size, n);
    ssize_t i;
    
    for (i = 0; i < nsize; i++) {
        array_set(a, i, val);
    }
    
    for(; i < size; i++) {
        array_set(a, i, val0);
    }
    
    return (void *)val + a->elt_size;
}


void * array_assign_array (struct array *a, const void *ptr, ssize_t n)
{
    assert(a);
    assert(n >= 0);
    assert(ptr || n == 0);

    char val0[a->elt_size];
    
    if (n < a->size)
        memset(val0, 0, a->elt_size);
    
    return array_assign_array_with(a, ptr, n, val0);
}


void * array_assign_array_with (struct array *a, const void *ptr, ssize_t n,
                               const void *val0)
{
    assert(a);
    assert(n >= 0);
    assert(ptr || n == 0);
    assert(val0 || n >= array_size(a));
    
    ssize_t size = a->size;
    ssize_t nsize = MIN(size, n);
    ssize_t nbytes = nsize * a->elt_size;
    ssize_t i;
    
    memcpy(array_begin(a), ptr, nbytes);
    
    for (i = nsize; i < size; i++) {
        array_set(a, i, val0);
    }
    
    return (void *)ptr + nbytes;
}


void array_copy (const struct array *a, struct array *dst)
{
    assert(a);
    assert(dst);
    
    array_assign_array(dst, array_begin(a), array_size(a));
}

void array_copy_with (const struct array *a, struct array *dst,
                      const void *val0)
{
    assert(a);
    assert(dst);
    assert(val0 || array_size(a) >= array_size(dst));
    
    array_assign_array_with(dst, array_begin(a), array_size(a), val0);
}


void array_swap (struct array *a, struct array *b)
{
    assert(a);
    assert(b);
    assert(array_size(a) == array_size(b));
    
    ssize_t i, n = array_size(a);
    char tmp[a->elt_size];
    
    for (i = 0; i < n; i++) {
        array_get(a, i, tmp);             // tmp  := a[i]
        array_set(a, i, array_ptr(b, i)); // a[i] := b[i]
        array_set(b, i, tmp);             // b[i] := tmp
    }
}


ssize_t array_lfind (const struct array *a, const void *key,
                    compare_fn compar)
{
    assert(a);
    assert(key);
    assert(compar);
    
    char *ptr = array_begin(a);
    size_t elt_size = _array_elt_size(a);
    ssize_t n = array_size(a);
    ssize_t i;
    
    for (i = 0; i < n; ptr += elt_size, i++) {
        if (compar(key, ptr) == 0)
            return i;
    }
    
    return ~i;
}

ssize_t array_bsearch (const struct array *a, const void *key,
                       compare_fn compar)
{
    assert(a);
    assert(key);
    assert(compar);
    
    size_t elt_size = _array_elt_size(a);
    void   *base = array_begin(a);
    ssize_t begin = 0;
    ssize_t end = array_size(a);
    ssize_t i;
    void *ptr;
    int cmp;
    
    while (begin < end) {
        i = begin + ((end - begin) >> 1);
        ptr = base + i * elt_size;
        cmp = compar(key, ptr);
        
        if (cmp < 0) {         // value < array[i]
            end = i;
        } else if (cmp > 0) {  // value > array[i]
            begin = i + 1;
        } else {               // value == array[i]
            return i;
        }
    }
    
    return ~end;               // begin == end, not found
}
