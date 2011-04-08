#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "pqueue.h"

struct pqueue * _pqueue_init (struct pqueue *q, compare_fn compar,
                              size_t elt_size)
{
    assert(compar);    
    assert(elt_size > 0);

    
    if (_darray_init(&q->array, elt_size)) {
        q->compare = compar;
        return q;
    }
    
    return NULL;
}


struct pqueue * pqueue_init_copy (struct pqueue *q, const struct pqueue *src)
{
    assert(q);
    assert(src);
    
    if (darray_init_copy(&q->array, &src->array)) {
        q->compare = src->compare;
        return q;
    }

    return NULL;
}


void pqueue_deinit (struct pqueue *q)
{
    assert(q);
    darray_deinit(&q->array);
}


void pqueue_assign_copy (struct pqueue *q, const struct pqueue *src)
{
    assert(q);
    assert(src);
    
    darray_assign_copy(&q->array, &src->array);
    q->compare = src->compare;
}


void * pqueue_copy_to (const struct pqueue *q, void *dst)
{
    assert(q);
    assert(dst);
    return darray_copy_to(&q->array, dst);
}


void * pqueue_push (struct pqueue *q, const void *val)
{
    assert(q);
    assert(val);
    
    struct darray *array = &q->array;
    compare_fn compare = q->compare;
    ssize_t icur = darray_size(array);

    // make space for the new element;
    darray_resize(array, icur + 1);
    
    // while current element has a parent:
    while (icur > 0) {
        ssize_t iparent = (icur - 1) >> 1;
        void *parent = darray_ptr(array, iparent);
        
        // if cur <= parent, heap condition is satisfied
        if (compare(val, parent) <= 0)
            break;
        
        // otherwise, swap(cur,parent)
        darray_set(array, icur, parent);
        icur = iparent;
    }
    
    // actually copy new element
    return darray_set(array, icur, val);
}


void * pqueue_push_array (struct pqueue *q, const void *src, ssize_t n)
{
    assert(q);
    assert(src || n == 0);
    assert(n >= 0);
    
    ssize_t i;
    
    for (i = 0; i < n; i++) {
        src = pqueue_push(q, src);
    }

    return (void *)src;
}


void pqueue_pop (struct pqueue *q)
{
    assert(q);
    assert(!pqueue_empty(q));
    
    struct darray *array = &q->array;
    ssize_t n = darray_size(array) - 1;
    
    if (n == 0)
        goto out;
    
    // swap the last element in the tree with the root, then heapify
    compare_fn compare = q->compare;
    void *cur = darray_ptr(array, n);
    ssize_t icur = 0;
    
    // while current element has at least one child
    while ((icur << 1) + 1 < n) {
        ssize_t ileft = (icur << 1) + 1;
        ssize_t iright = (icur << 1) + 2;
        ssize_t imax;
        void *left = darray_ptr(array, ileft);
        void *right = darray_ptr(array, iright);
        void *max;

        // find the child with highest priority
        if (iright == n || compare(right, left) <= 0) {
            imax = ileft;
            max = left;
        } else {
            imax = iright;
            max = right;
        }
        
        // stop if heap condition is satisfied
        if (compare(max, cur) <= 0)
            break;
            
        // otherwise swap current with maximum child
        darray_set(array, icur, max);
        icur = imax;
    }
    
    // actually do the copy
    darray_set(array, icur, cur);
    
out:
    darray_resize(array, n);
}


void * pqueue_pop_array (struct pqueue *q, void *dst, ssize_t n)
{
    assert(q);
    assert(dst || n == 0);
    assert(n >= 0);
    assert(n <= pqueue_size(q));
    
    size_t i;
    
    for (i = 0; i < n; i++) {
        dst = pqueue_get_top(q, dst);
        pqueue_pop(q);
    }

    return dst;
}
