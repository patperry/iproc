#ifndef _PQUEUE_H
#define _PQUEUE_H

/* Priority Queue type.
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include "array.h"
#include "compare.h"

struct pqueue {
	struct array array;
	compare_fn compare;
};

/* create, destroy */
void pqueue_init(struct pqueue *q, compare_fn compar, size_t elt_size);
void pqueue_init_copy(struct pqueue *q, const struct pqueue *src);
void pqueue_assign_copy(struct pqueue *q, const struct pqueue *src);
void pqueue_deinit(struct pqueue *q);


void pqueue_copy_to(const struct pqueue *q, void *dst);
void pqueue_clear(struct pqueue *q);

/* informative */
static inline void *pqueue_top(const struct pqueue *q);
static inline ssize_t pqueue_count(const struct pqueue *q);
static inline size_t pqueue_elt_size(const struct pqueue *q);

/* operations */
bool pqueue_push(struct pqueue *q, const void *val);
void pqueue_pop(struct pqueue *q);
void pqueue_update_top(struct pqueue *q);
bool pqueue_set_capacity(struct pqueue *q, ssize_t n);

/* inline function definitions */
ssize_t pqueue_count(const struct pqueue *q)
{
	return array_count(&q->array);
}

size_t pqueue_elt_size(const struct pqueue *q)
{
	return array_elt_size(&q->array);
}

void *pqueue_top(const struct pqueue *q)
{
	assert(pqueue_count(q));
	return array_item(&q->array, 0);
}

#endif /* _PQUEUE_H */
