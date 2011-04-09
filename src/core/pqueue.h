#ifndef _PQUEUE_H
#define _PQUEUE_H

/* Priority Queue type.
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include "darray.h"
#include "compare.h"


struct pqueue {
    struct darray  array;
    compare_fn     compare;
};


/* create, destroy */
#define         pqueue_init(q,t,compar) _pqueue_init(q,compar,sizeof(t))
struct pqueue * pqueue_init_copy (struct pqueue *q, const struct pqueue *src);
void            pqueue_deinit    (struct pqueue *q);

/* assignment, copy */
struct pqueue * pqueue_assign_copy (struct pqueue *q, const struct pqueue *src);
void * pqueue_copy_to     (const struct pqueue *q, void *dst);

/* informative */
#define                    pqueue_top(q,t) (*((t const * const)_pqueue_top(q)))
static inline void         pqueue_get_top  (const struct pqueue *q, void *dst);
static inline bool         pqueue_empty    (const struct pqueue *q);
static inline ssize_t      pqueue_size     (const struct pqueue *q);
static inline ssize_t      pqueue_max_size (const struct pqueue *q);
static inline size_t       pqueue_elt_size (const struct pqueue *q);

/* operations */
void   pqueue_push       (struct pqueue *q, const void *val);
void * pqueue_push_array (struct pqueue *q, const void *src, ssize_t n);
void   pqueue_pop        (struct pqueue *q);
void * pqueue_pop_array  (struct pqueue *q, void *dst, ssize_t n);



/* private functions */
struct pqueue * _pqueue_init (struct pqueue *q, compare_fn compar,
                              size_t elt_size);
static inline const void * _pqueue_top (const struct pqueue *q);


/* inline function definitions */
bool    pqueue_empty    (const struct pqueue *q) { return darray_empty(&q->array); }
ssize_t pqueue_size     (const struct pqueue *q) { return darray_size(&q->array); }
ssize_t pqueue_max_size (const struct pqueue *q) { return darray_max_size(&q->array); }
size_t  pqueue_elt_size (const struct pqueue *q) { return darray_elt_size(&q->array); }


const void * _pqueue_top (const struct pqueue *q)
{
    assert(!pqueue_empty(q));
    return darray_begin(&q->array);
}

void pqueue_get_top (const struct pqueue *q, void *dst)
{
    assert(!pqueue_empty(q));
    return darray_get(&q->array, 0, dst);
}


#endif /* _PQUEUE_H */