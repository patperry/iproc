#include "port.h"

#include <assert.h>
#include <string.h>

#include "memory.h"
#include "pqueue.h"

static void
iproc_pqueue_free (iproc_pqueue *pqueue)
{
    if (pqueue) {
        darray_free(pqueue->array);
        iproc_free(pqueue);
    }
}

iproc_pqueue *
iproc_pqueue_new (size_t            eltsize,
                  iproc_compare_fn  compare)
{
    assert(eltsize > 0);
    assert(compare);
    
    iproc_pqueue *pqueue = iproc_malloc(sizeof(*pqueue));
    if (!pqueue)
        return NULL;
    
    pqueue->array = _darray_new(eltsize);
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
    
    result->array = darray_new_copy(pqueue->array);
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
    return darray_empty(pqueue->array);
}

ssize_t
iproc_pqueue_size (iproc_pqueue *pqueue)
{
    assert(pqueue);
    return darray_size(pqueue->array);
}


void *
iproc_pqueue_top (iproc_pqueue *pqueue)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    
    return darray_begin(pqueue->array);
}

void
iproc_pqueue_push (iproc_pqueue *pqueue,
                   void         *eltp)
{
    assert(pqueue);
    assert(eltp);
    
    struct darray *array = pqueue->array;
    size_t elt_size = darray_elt_size(array);
    iproc_compare_fn compare = pqueue->compare;
    ssize_t icur = darray_size(array);

    // make space for the new element;
    darray_resize(array, icur + 1);
    
    // while current element has a parent:
    while (icur > 0) {
        ssize_t iparent = (icur - 1) >> 1;
        void *parent = &darray_index(array, char, iparent * elt_size);
        
        // if cur <= parent, heap condition is satisfied
        if (compare(eltp, parent) <= 0)
            break;
        
        // otherwise, swap(cur,parent)
        darray_set(array, icur, parent);
        icur = iparent;
    }
    
    // actually copy new element
    darray_set(array, icur, eltp);
}

void
iproc_pqueue_push_array (iproc_pqueue *pqueue,
                         void         *elts,
                         ssize_t       n)
{
    assert(pqueue);
    assert(elts || n == 0);
    assert(n >= 0);
    
    size_t elem_size = darray_elt_size(pqueue->array);
    void *end = elts + n * elem_size;
    
    for (; elts < end; elts += elem_size) {
        iproc_pqueue_push(pqueue, elts);
    }
}


void
iproc_pqueue_pop (iproc_pqueue *pqueue,
                  void         *eltp)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    assert(eltp);
    
    struct darray *array = pqueue->array;
    ssize_t n = darray_size(array) - 1;
    size_t elt_size = darray_elt_size(pqueue->array);
    
    // copy the top element
    memcpy(eltp, iproc_pqueue_top(pqueue), elt_size);
    
    if (n == 0)
        goto resize;
    
    // swap the last element in the tree with the root, then heapify
    iproc_compare_fn compare = pqueue->compare;
    void *cur = &darray_index(array, char, n * elt_size);
    ssize_t icur = 0;
    
    // while current element has at least one child
    while ((icur << 1) + 1 < n) {
        ssize_t ileft = (icur << 1) + 1;
        ssize_t iright = (icur << 1) + 2;
        ssize_t imax;
        void *left = &darray_index(array, char, ileft * elt_size);
        void *right = &darray_index(array, char, iright * elt_size);
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
    
resize:
    darray_resize(array, n);
}

void
iproc_pqueue_pop_array (iproc_pqueue *pqueue,
                        void         *elts,
                        ssize_t       n)
{
    assert(pqueue);
    assert(elts || n == 0);
    assert(n >= 0);
    assert(n <= iproc_pqueue_size(pqueue));
    
    size_t elt_size = darray_elt_size(pqueue->array);
    void *end = elts + n * elt_size;
    
    for (; elts < end; elts += elt_size) {
        iproc_pqueue_pop(pqueue, elts);
    }
}
