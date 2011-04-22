#ifndef _ARRAY_H
#define _ARRAY_H

/* Generic Array type.
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include <string.h>		// memcpy
#include "compare.h"

struct array {
	void *data;
	ssize_t size;
	size_t elt_size;
	bool owner;
};

/* create, destroy */
bool array_init(struct array *a, ssize_t size, size_t elt_size);
bool array_init_copy(struct array *a, const struct array *src);
bool array_reinit(struct array *a, ssize_t n, size_t elt_size);
void array_deinit(struct array *a);

/* views */
void array_init_view(struct array *a,
		     const void *data, ssize_t size, size_t elt_size);
void array_init_slice(struct array *a,
		      const struct array *parent, ssize_t i, ssize_t n);

/* assign, copy, fill */
void array_assign_array(struct array *a, const void *src);
void array_assign_copy(struct array *a, const struct array *src);
void *array_copy_to(const struct array *a, void *dst);
void *array_copy_range_to(const struct array *a, ssize_t i, ssize_t n, void *dst);
void array_fill(struct array *a, const void *val);
void array_fill_range(struct array *a, ssize_t i, ssize_t n, const void *val);

/* index */
static inline void *array_at(const struct array *a, ssize_t i);
static inline void array_set(struct array *a, ssize_t i, const void *src);
static inline void array_set_range(struct array *a,
				   ssize_t i, ssize_t n, const void *src);
#define array_front(a)           (array_at(a, 0))
#define array_set_front(a, val)  (array_set(a, 0, val))
#define array_back(a)            (array_at(a, array_size(a) - 1))
#define array_set_back(a, val)   (array_set(a, array_size(a) - 1, val))


/* informative */
static inline ssize_t array_size(const struct array *a);
static inline bool array_empty(const struct array *a);
static inline size_t array_elt_size(const struct array *a);
static inline ssize_t array_max_size(const struct array *a);
static inline bool array_owner(const struct array *a);
static inline bool array_overlaps(const struct array *a,
				  ssize_t i,
				  ssize_t n, const void *ptr, ssize_t nel);

/* operations */
void array_swap(struct array *a, ssize_t i, ssize_t j);
void array_reverse(struct array *a);

/* searching, sorting */
bool array_contains(const struct array *a, const void *key, equals_fn equal);
void *array_find(const struct array *a, const void *key, equals_fn equal);
ssize_t array_find_index(const struct array *a,
			 const void *key, equals_fn equal);
void *array_find_last(const struct array *a, const void *key, equals_fn equal);
ssize_t array_find_last_index(const struct array *a,
			      const void *key, equals_fn equal);
ssize_t array_binary_search(const struct array *a,
			    const void *key, compare_fn compar);
void array_sort(struct array *a, compare_fn compar);


/* inline function defs */
ssize_t array_size(const struct array *a)
{
	return a->size;
}

bool array_empty(const struct array * a)
{
	return a->size == 0;
}

size_t array_elt_size(const struct array * a)
{
	return a->elt_size;
}

ssize_t array_max_size(const struct array * a)
{
	return SSIZE_MAX / a->elt_size;
}

bool array_owner(const struct array * a)
{
	return a->owner;
}

bool array_overlaps(const struct array * a, ssize_t i, ssize_t n,
		    const void *ptr, ssize_t nel)
{
	assert(a);
	assert(n >= 0);
	assert(0 <= i && i <= array_size(a) - n);
	assert(0 <= nel && nel < array_max_size(a));

	if (array_empty(a) || n == 0 || nel == 0)
		return false;
	
	size_t elt_size = array_elt_size(a);
	const void *begin1 = array_at(a, i);
	const void *end1 = begin1 + n * elt_size;
	const void *begin2 = ptr;
	const void *end2 = begin2 + nel * elt_size;

	return ((begin1 <= begin2 && begin2 < end1)
		|| (begin2 <= begin1 && begin1 < end2));
}

void *array_at(const struct array *a, ssize_t i)
{
	assert(0 <= i && i < array_size(a));
	return a->data + i * a->elt_size;
	
}

void array_set(struct array *a, ssize_t i, const void *src)
{
	assert(0 <= i && i < array_size(a));
	assert(!array_overlaps(a, i, 1, src, 1));

	array_set_range(a, i, 1, src);
}

void array_set_range(struct array *a, ssize_t i, ssize_t n, const void *src)
{
	assert(n >= 0);
	assert(0 <= i && i <= array_size(a) - n);
	assert(!array_overlaps(a, i, n, src, n));

	if (n == 0)
		return;
	
	size_t nbytes = n * array_elt_size(a);
	memcpy(array_at(a, i), src, nbytes);
}

#endif /* _ARRAY_H */
