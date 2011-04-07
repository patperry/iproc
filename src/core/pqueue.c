#include "port.h"

#include <assert.h>
#include <string.h>

#include "memory.h"
#include "pqueue.h"

static void
iproc_pqueue_free (iproc_pqueue *pqueue)
{
    if (pqueue) {
        darray_deinit(&pqueue->array);
        iproc_free(pqueue);
    }
}

iproc_pqueue *
iproc_pqueue_new (size_t     eltsize,
                  compare_fn compare)
{
    assert(eltsize > 0);
    assert(compare);
    
    iproc_pqueue *pqueue = iproc_malloc(sizeof(*pqueue));
    if (!pqueue)
        return NULL;
    
    _darray_init(&pqueue->array, eltsize);
    pqueue->compare = compare;
    iproc_refcount_init(&pqueue->refcount);
    
    return pqueue;
}

iproc_pqueue *
iproc_pqueue_new_copy (iproc_pqueue *pqueue)
{
    assert(pqueue);
    
    iproc_pqueue *result = iproc_malloc(sizeof(*result));
    if (!result)
        return NULL;
    
    darray_init_copy(&result->array, &pqueue->array);
    result->compare = pqueue->compare;
    iproc_refcount_init(&result->refcount);

    return result;
}

void
iproc_pqueue_ref (iproc_pqueue *pqueue)
{
    if (pqueue) {
        iproc_refcount_get(&pqueue->refcount);
    }
}


static void
iproc_pqueue_release (iproc_refcount *refcount)
{
    iproc_pqueue *pqueue = container_of(refcount, iproc_pqueue, refcount);
    iproc_pqueue_free(pqueue);
}

void
iproc_pqueue_unref (iproc_pqueue *pqueue)
{
    if (!pqueue)
        return;
    
    iproc_refcount_put(&pqueue->refcount, iproc_pqueue_release);
}

bool
iproc_pqueue_empty (iproc_pqueue *pqueue)
{
    assert(pqueue);
    return darray_empty(&pqueue->array);
}

ssize_t
iproc_pqueue_size (iproc_pqueue *pqueue)
{
    assert(pqueue);
    return darray_size(&pqueue->array);
}


const void *
iproc_pqueue_top (iproc_pqueue *pqueue)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    
    return darray_begin(&pqueue->array);
}

void *
iproc_pqueue_get_top (iproc_pqueue *pqueue, void *dst)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    assert(dst);
    
    return darray_get(&pqueue->array, 0, dst);
}

void *
iproc_pqueue_push (iproc_pqueue *pqueue,
                   const void   *eltp)
{
    assert(pqueue);
    assert(eltp);
    
    struct darray *array = &pqueue->array;
    compare_fn compare = pqueue->compare;
    ssize_t icur = darray_size(array);

    // make space for the new element;
    darray_resize(array, icur + 1);
    
    // while current element has a parent:
    while (icur > 0) {
        ssize_t iparent = (icur - 1) >> 1;
        void *parent = darray_ptr(array, iparent);
        
        // if cur <= parent, heap condition is satisfied
        if (compare(eltp, parent) <= 0)
            break;
        
        // otherwise, swap(cur,parent)
        darray_set(array, icur, parent);
        icur = iparent;
    }
    
    // actually copy new element
    return darray_set(array, icur, eltp);
}

void *
iproc_pqueue_push_array (iproc_pqueue *pqueue,
                         const void   *elts,
                         ssize_t       n)
{
    assert(pqueue);
    assert(elts || n == 0);
    assert(n >= 0);
    
    ssize_t i;
    
    for (i = 0; i < n; i++) {
        elts = iproc_pqueue_push(pqueue, elts);
    }

    return (void *)elts;
}


void
iproc_pqueue_pop (iproc_pqueue *pqueue)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    
    struct darray *array = &pqueue->array;
    ssize_t n = darray_size(array) - 1;
    
    if (n == 0)
        goto out;
    
    // swap the last element in the tree with the root, then heapify
    compare_fn compare = pqueue->compare;
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

void *
iproc_pqueue_pop_array (iproc_pqueue *pqueue,
                        void         *elts,
                        ssize_t       n)
{
    assert(pqueue);
    assert(elts || n == 0);
    assert(n >= 0);
    assert(n <= iproc_pqueue_size(pqueue));
    
    size_t i;
    
    for (i = 0; i < n; i++) {
        elts = iproc_pqueue_get_top(pqueue, elts);
        iproc_pqueue_pop(pqueue);
    }

    return elts;
}
