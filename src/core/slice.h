#ifndef _SLICE_H
#define _SLICE_H

/* Generic Slice type.
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include <string.h>		// memcpy
#include "compare.h"
#include "delegate.h"
#include "util.h"

struct slice {
	void *data;
	ssize_t size;
	size_t elt_size;
};

static inline struct slice slice_make(const void *ptr, ssize_t n,
				      size_t elt_size);
static inline void slice_init(struct slice *a, const void *data, ssize_t size,
			      size_t elt_size);

/* assign, copy, fill */
void slice_copy_to(const struct slice *a, void *dst);
void slice_fill(struct slice *a, const void *val);

/* properties */
static inline ssize_t slice_count(const struct slice *a);
static inline void *slice_item(const struct slice *a, ssize_t i);
static inline void slice_set_item(struct slice *a, ssize_t i, const void *src);
static inline size_t slice_elt_size(const struct slice *a);

/* operations */
void slice_reverse(struct slice *a);

/* searching, sorting */
bool slice_exists(const struct slice *a, predicate_fn match, void *udata);
void *slice_find(const struct slice *a, predicate_fn match, void *udata);
ssize_t slice_find_index(const struct slice *a, predicate_fn match,
			 void *udata);
void *slice_find_last(const struct slice *a, predicate_fn match, void *udata);
ssize_t slice_find_last_index(const struct slice *a, predicate_fn match,
			      void *udata);
ssize_t slice_binary_search(const struct slice *a, const void *key,
			    compare_fn compar);
void slice_sort(struct slice *a, compare_fn compar);

/* inline function defs */
struct slice slice_make(const void *ptr, ssize_t n, size_t elt_size)
{
	assert(ptr || n == 0);
	assert(0 <= n && n <= SSIZE_MAX / MAX(1, elt_size));
	assert(elt_size >= 0);

	struct slice a;
	slice_init(&a, ptr, n, elt_size);
	return a;
}

void slice_init(struct slice *a, const void *ptr, ssize_t n, size_t elt_size)
{
	assert(a);
	assert(ptr || n == 0);
	assert(0 <= n && n <= SSIZE_MAX / MAX(1, elt_size));
	assert(elt_size >= 0);

	a->data = (void *)ptr;
	a->size = n;
	a->elt_size = elt_size;
}

ssize_t slice_count(const struct slice *a)
{
	assert(a);
	return a->size;
}

size_t slice_elt_size(const struct slice *a)
{
	assert(a);
	return a->elt_size;
}

void *slice_item(const struct slice *a, ssize_t i)
{
	assert(a);
	assert(0 <= i && i < slice_count(a));
	return (char *)a->data + i * a->elt_size;

}

void slice_set_item(struct slice *a, ssize_t i, const void *src)
{
	assert(a);
	assert(0 <= i && i < slice_count(a));
	assert(!memory_overlaps
	       (slice_item(a, i), 1, src, 1, slice_elt_size(a)));

	memcpy(slice_item(a, i), src, slice_elt_size(a));
}

#endif /* _SLICE_H */
