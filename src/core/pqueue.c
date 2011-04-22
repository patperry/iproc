#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "pqueue.h"

bool pqueue_init(struct pqueue *q, compare_fn compar, size_t elt_size)
{
	assert(q);
	assert(compar);
	assert(elt_size > 0);

	if (darray_init(&q->array, elt_size)) {
		q->compare = compar;
		return q;
	}

	return NULL;
}

bool pqueue_init_copy(struct pqueue *q, const struct pqueue *src)
{
	assert(q);
	assert(src);

	if (darray_init_copy(&q->array, &src->array)) {
		q->compare = src->compare;
		return q;
	}

	return NULL;
}

void pqueue_deinit(struct pqueue *q)
{
	assert(q);
	darray_deinit(&q->array);
}

struct pqueue *pqueue_assign_copy(struct pqueue *q, const struct pqueue *src)
{
	assert(q);
	assert(src);

	if (darray_assign_copy(&q->array, &src->array)) {
		q->compare = src->compare;
		return q;
	}

	return NULL;
}

void *pqueue_copy_to(const struct pqueue *q, void *dst)
{
	assert(q);
	assert(dst);
	return darray_copy_to(&q->array, dst);
}

void *pqueue_push(struct pqueue *q, const void *val)
{
	assert(q);
	assert(val);

	struct darray *array = &q->array;
	compare_fn compare = q->compare;
	ssize_t icur = darray_size(array);

	// make space for the new element
	if (!darray_reserve(array, icur + 1))
		return NULL;

	darray_resize(array, icur + 1);

	// while current element has a parent:
	while (icur > 0) {
		ssize_t iparent = (icur - 1) >> 1;
		void *parent = darray_at(array, iparent);

		// if cur <= parent, heap condition is satisfied
		if (compare(val, parent) <= 0)
			break;

		// otherwise, swap(cur,parent)
		darray_set(array, icur, parent);
		icur = iparent;
	}

	// actually copy new element
	darray_set(array, icur, val);

	return (void *)val + pqueue_elt_size(q);
}

void *pqueue_push_all(struct pqueue *q, const void *src, ssize_t n)
{
	assert(q);
	assert(src || n == 0);
	assert(n >= 0);

	// make space for the new elements
	if (!darray_reserve(&q->array, pqueue_size(q) + n))
		return NULL;

	size_t elt_size = pqueue_elt_size(q);
	const void *end = src + n * elt_size;

	for (; src < end; src += elt_size) {
		pqueue_push(q, src);
	}

	return (void *)src;
}

void pqueue_pop(struct pqueue *q)
{
	assert(q);
	assert(!pqueue_empty(q));

	struct darray *array = &q->array;
	ssize_t n = darray_size(array) - 1;

	if (n == 0)
		goto out;

	// swap the last element in the tree with the root, then heapify
	compare_fn compare = q->compare;
	void *cur = darray_back(array);
	ssize_t icur = 0;

	// while current element has at least one child
	while ((icur << 1) + 1 < n) {
		ssize_t ileft = (icur << 1) + 1;
		ssize_t iright = (icur << 1) + 2;
		ssize_t imax;
		void *left = darray_at(array, ileft);
		void *right = darray_at(array, iright);
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
