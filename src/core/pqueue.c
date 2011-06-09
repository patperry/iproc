#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "pqueue.h"

void pqueue_init(struct pqueue *q, compare_fn compar, size_t elt_size)
{
	assert(q);
	assert(compar);
	assert(elt_size > 0);

	array_init(&q->array, elt_size);
	q->compare = compar;
}

void pqueue_init_copy(struct pqueue *q, const struct pqueue *src)
{
	assert(q);
	assert(src);

	array_init_copy(&q->array, &src->array);
	q->compare = src->compare;
}

void pqueue_deinit(struct pqueue *q)
{
	assert(q);
	array_deinit(&q->array);
}

void pqueue_assign_copy(struct pqueue *q, const struct pqueue *src)
{
	assert(q);
	assert(src);

	array_assign_copy(&q->array, &src->array);
	q->compare = src->compare;
}

void pqueue_copy_to(const struct pqueue *q, void *dst)
{
	assert(q);
	assert(dst);
	array_copy_to(&q->array, dst);
}

void pqueue_clear(struct pqueue *q)
{
	assert(q);
	array_clear(&q->array);
}

bool pqueue_push(struct pqueue *q, const void *val)
{
	assert(q);
	assert(val);

	struct array *array = &q->array;
	compare_fn compare = q->compare;
	ssize_t icur = array_count(array);

	// make space for the new element
	array_add(array, NULL);

	// while current element has a parent:
	while (icur > 0) {
		ssize_t iparent = (icur - 1) >> 1;
		void *parent = array_item(array, iparent);

		// if cur <= parent, heap condition is satisfied
		if (compare(val, parent) <= 0)
			break;

		// otherwise, swap(cur,parent)
		array_set_item(array, icur, parent);
		icur = iparent;
	}

	// actually copy new element
	array_set_item(array, icur, val);

	return true;
}

void pqueue_pop(struct pqueue *q)
{
	assert(q);
	assert(pqueue_count(q));

	struct array *array = &q->array;
	ssize_t n = array_count(array) - 1;

	if (n == 0)
		goto out;

	// swap the last element in the tree with the root, then heapify
	compare_fn compare = q->compare;
	void *cur = array_item(array, array_count(array) - 1);
	ssize_t icur = 0;

	// while current element has at least one child
	while ((icur << 1) + 1 < n) {
		ssize_t ileft = (icur << 1) + 1;
		ssize_t iright = (icur << 1) + 2;
		ssize_t imax;
		void *left = array_item(array, ileft);
		void *right = array_item(array, iright);
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
		array_set_item(array, icur, max);
		icur = imax;
	}

	// actually do the copy
	array_set_item(array, icur, cur);

out:
	array_remove_at(array, n);
}

void pqueue_update_top(struct pqueue *q)
{
	assert(q);
	assert(pqueue_count(q));

	struct array *array = &q->array;
	ssize_t n = array_count(array);
	size_t elt_size = pqueue_elt_size(q);
	compare_fn compare = q->compare;

	// make a temporary copy of the old top
	void *cur = alloca(elt_size);
	memcpy(cur, array_item(array, 0), elt_size);
	ssize_t icur = 0;

	// while current element has at least one child
	while ((icur << 1) + 1 < n) {
		ssize_t ileft = (icur << 1) + 1;
		ssize_t imax;
		void *left = array_item(array, ileft);
		void *max;

		// find the child with highest priority
		if (ileft == n - 1) {
			imax = ileft;
			max = left;
		} else {
			ssize_t iright = ileft + 1;
			void *right = array_item(array, iright);

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
		array_set_item(array, icur, max);
		icur = imax;
	}

	// actually do the copy
	array_set_item(array, icur, cur);
}

bool pqueue_set_capacity(struct pqueue *q, ssize_t n)
{
	assert(q);
	assert(n >= 0);
	array_set_capacity(&q->array, n);
	return true;
}
