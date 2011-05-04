#ifndef _PQUEUE_H
#define _PQUEUE_H

/* Priority Queue type.
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include "darray.h"
#include "compare.h"

struct pqueue {
	struct darray array;
	compare_fn compare;
};

/* create, destroy */
bool pqueue_init(struct pqueue *q, compare_fn compar, size_t elt_size);
bool pqueue_init_copy(struct pqueue *q, const struct pqueue *src);
void pqueue_deinit(struct pqueue *q);

/* assignment, copy, clear */
struct pqueue *pqueue_assign_copy(struct pqueue *q, const struct pqueue *src);
void *pqueue_copy_to(const struct pqueue *q, void *dst);
void pqueue_clear(struct pqueue *q);

/* informative */
static inline void *pqueue_top(const struct pqueue *q);
static inline bool pqueue_empty(const struct pqueue *q);
static inline ssize_t pqueue_size(const struct pqueue *q);
static inline ssize_t pqueue_max_size(const struct pqueue *q);
static inline size_t pqueue_elt_size(const struct pqueue *q);

/* operations */
bool pqueue_push(struct pqueue *q, const void *val);
bool pqueue_push_all(struct pqueue *q, const void *vals, ssize_t n);
void pqueue_pop(struct pqueue *q);
void pqueue_update_top(struct pqueue *q);
bool pqueue_reserve(struct pqueue *q, ssize_t n);
bool pqueue_reserve_push(struct pqueue *q, ssize_t npush);


/* inline function definitions */
bool pqueue_empty(const struct pqueue *q)
{
	return darray_empty(&q->array);
}

ssize_t pqueue_size(const struct pqueue *q)
{
	return darray_size(&q->array);
}

ssize_t pqueue_max_size(const struct pqueue *q)
{
	return darray_max_size(&q->array);
}

size_t pqueue_elt_size(const struct pqueue *q)
{
	return darray_elt_size(&q->array);
}

void *pqueue_top(const struct pqueue *q)
{
	assert(!pqueue_empty(q));
	return darray_front(&q->array);
}

#endif /* _PQUEUE_H */
