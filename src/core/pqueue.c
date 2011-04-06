#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <assert.h>

#include "memory.h"
#include "pqueue.h"

static void
iproc_pqueue_free (iproc_pqueue *pqueue)
{
    if (pqueue) {
        iproc_array_unref(pqueue->array);
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
    
    pqueue->array = iproc_array_new(eltsize);
    pqueue->compare = compare;
    iproc_refcount_init(&pqueue->refcount);
    
    return pqueue;
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
    return iproc_array_empty(pqueue->array);
}

ssize_t
iproc_pqueue_size (iproc_pqueue *pqueue)
{
    assert(pqueue);
    return iproc_array_size(pqueue->array);
}


void *
iproc_pqueue_top (iproc_pqueue *pqueue)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    
    return &iproc_array_index(pqueue->array, char, 0);
}

void
iproc_pqueue_push (iproc_pqueue *pqueue,
                   void         *eltp)
{
    assert(pqueue);
    assert(eltp);
    
    iproc_array *array = pqueue->array;
    size_t eltsize = iproc_array_elem_size(array);
    iproc_compare_fn compare = pqueue->compare;
    ssize_t icur = iproc_array_size(array);

    // make space for the new element;
    iproc_array_set_size(array, icur + 1);
    
    // while current element has a parent:
    while (icur > 0) {
        ssize_t iparent = (icur - 1) >> 1;
        void *parent = &iproc_array_index(array, char, iparent * eltsize);
        
        // if cur <= parent, heap condition is satisfied
        if (compare(eltp, parent) <= 0)
            break;
        
        // otherwise, swap(cur,parent)
        iproc_array_set(array, icur, parent);
        icur = iparent;
    }
    
    // actually copy new element
    iproc_array_set(array, icur, eltp);
}

void
iproc_pqueue_push_array (iproc_pqueue *pqueue,
                         void         *elts,
                         ssize_t       n)
{
    assert(pqueue);
    assert(n == 0 || elts);
    assert(n >= 0);
    
    size_t elem_size = pqueue->array->elem_size;
    void *end = elts + n * elem_size;
    
    for (; elts < end; elts += elem_size) {
        iproc_pqueue_push(pqueue, elts);
    }
}


void
iproc_pqueue_pop (iproc_pqueue *pqueue)
{
    assert(pqueue);
    assert(!iproc_pqueue_empty(pqueue));
    
    iproc_array *array = pqueue->array;
    ssize_t n = iproc_array_size(array) - 1;
    
    if (n == 0)
        goto resize;
    
    // swap the last element in the tree with the root, then heapify
    size_t eltsize = iproc_array_elem_size(array);
    iproc_compare_fn compare = pqueue->compare;
    void *cur = &iproc_array_index(array, char, n * eltsize);
    ssize_t icur = 0;
    
    // while current element has at least one child
    while ((icur << 1) + 1 < n) {
        ssize_t ileft = (icur << 1) + 1;
        ssize_t iright = (icur << 1) + 2;
        ssize_t imax;
        void *left = &iproc_array_index(array, char, ileft * eltsize);
        void *right = &iproc_array_index(array, char, iright * eltsize);
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
        iproc_array_set(array, icur, max);
        icur = imax;
    }
    
    // actually do the copy
    iproc_array_set(array, icur, cur);
    
resize:
    iproc_array_set_size(array, n);
}
