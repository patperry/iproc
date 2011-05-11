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

	if (list_init(&q->array, elt_size)) {
		q->compare = compar;
		return q;
	}

	return NULL;
}

bool pqueue_init_copy(struct pqueue *q, const struct pqueue *src)
{
	assert(q);
	assert(src);

	if (list_init_copy(&q->array, &src->array)) {
		q->compare = src->compare;
		return q;
	}

	return NULL;
}

void pqueue_deinit(struct pqueue *q)
{
	assert(q);
	list_deinit(&q->array);
}

struct pqueue *pqueue_assign_copy(struct pqueue *q, const struct pqueue *src)
{
	assert(q);
	assert(src);

	if (list_assign_copy(&q->array, &src->array)) {
		q->compare = src->compare;
		return q;
	}

	return NULL;
}

void pqueue_copy_to(const struct pqueue *q, void *dst)
{
	assert(q);
	assert(dst);
	list_copy_to(&q->array, dst);
}

void pqueue_clear(struct pqueue *q)
{
	assert(q);
	list_clear(&q->array);
}

bool pqueue_push(struct pqueue *q, const void *val)
{
	assert(q);
	assert(val);

	struct list *array = &q->array;
	compare_fn compare = q->compare;
	ssize_t icur = list_count(array);

	// make space for the new element
	if (!list_set_capacity(array, icur + 1))
		return false;

	list_add(array, NULL);

	// while current element has a parent:
	while (icur > 0) {
		ssize_t iparent = (icur - 1) >> 1;
		void *parent = list_item(array, iparent);

		// if cur <= parent, heap condition is satisfied
		if (compare(val, parent) <= 0)
			break;

		// otherwise, swap(cur,parent)
		list_set_item(array, icur, parent);
		icur = iparent;
	}

	// actually copy new element
	list_set_item(array, icur, val);

	return true;
}

bool pqueue_push_all(struct pqueue *q, const void *src, ssize_t n)
{
	assert(q);
	assert(src || n == 0);
	assert(n >= 0);

	// make space for the new elements
	if (!pqueue_reserve_push(q, n))
		return false;

	size_t elt_size = pqueue_elt_size(q);
	const void *end = (char *)src + n * elt_size;

	for (; src < end; src = (char *)src + elt_size) {
		pqueue_push(q, src);
	}

	return true;
}

void pqueue_pop(struct pqueue *q)
{
	assert(q);
	assert(!pqueue_empty(q));

	struct list *array = &q->array;
	ssize_t n = list_count(array) - 1;

	if (n == 0)
		goto out;

	// swap the last element in the tree with the root, then heapify
	compare_fn compare = q->compare;
	void *cur = list_item(array, list_count(array) - 1);
	ssize_t icur = 0;

	// while current element has at least one child
	while ((icur << 1) + 1 < n) {
		ssize_t ileft = (icur << 1) + 1;
		ssize_t iright = (icur << 1) + 2;
		ssize_t imax;
		void *left = list_item(array, ileft);
		void *right = list_item(array, iright);
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
		list_set_item(array, icur, max);
		icur = imax;
	}

	// actually do the copy
	list_set_item(array, icur, cur);

out:
	list_remove_at(array, n);
}

void pqueue_update_top(struct pqueue *q)
{
	assert(q);
	assert(!pqueue_empty(q));

	struct list *array = &q->array;
	ssize_t n = list_count(array);
	size_t elt_size = pqueue_elt_size(q);
	compare_fn compare = q->compare;

	// make a temporary copy of the old top
	void *cur = alloca(elt_size);
	memcpy(cur, list_item(array, 0), elt_size);
	ssize_t icur = 0;

	// while current element has at least one child
	while ((icur << 1) + 1 < n) {
		ssize_t ileft = (icur << 1) + 1;
		ssize_t imax;
		void *left = list_item(array, ileft);
		void *max;

		// find the child with highest priority
		if (ileft == n - 1) {
			imax = ileft;
			max = left;
		} else {
			ssize_t iright = ileft + 1;
			void *right = list_item(array, iright);

			if (compare(right, left) <= 0) {
				imax = ileft;
				max = left;
			} else {
				imax = iright;
				max = right;
			}
		}

		// stop if heap condition is satisfied
		if (compare(max, cur) <= 0)
			break;

		// otherwise swap current with maximum child
		list_set_item(array, icur, max);
		icur = imax;
	}

	// actually do the copy
	list_set_item(array, icur, cur);
}

bool pqueue_reserve(struct pqueue *q, ssize_t n)
{
	assert(q);
	assert(n >= 0);
	return list_set_capacity(&q->array, n);
}

bool pqueue_reserve_push(struct pqueue *q, ssize_t npush)
{
	assert(q);
	assert(npush >= 0);

	return pqueue_reserve(q, pqueue_size(q) + npush);
}
