#ifndef _DARRAY_H
#define _DARRAY_H

/* Dynamic Array type
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include <string.h>		// memcpy
#include "array.h"

struct darray {
	struct array array;
	ssize_t size;
};

/* create, destroy */
struct darray *darray_init(struct darray *a, size_t elt_size);
bool darray_init_copy(struct darray *a, const struct darray *src);
void darray_deinit(struct darray *a);

/* assign, copy, fill */
bool darray_assign(struct darray *a, const void *ptr, ssize_t n);
bool darray_assign_repeat(struct darray *a, ssize_t n, const void *val);
bool darray_assign_copy(struct darray *a, const struct darray *src);
void *darray_copy_to(const struct darray *a, void *dst);
void *darray_copy_range_to(const struct darray *a, ssize_t i, ssize_t n,
			   void *dst);
void darray_fill(struct darray *a, const void *val);
void darray_fill_range(struct darray *a, ssize_t i, ssize_t n, const void *val);

/* index */
static inline void *darray_at(const struct darray *a, ssize_t i);
static inline void darray_set(struct darray *a, ssize_t i, const void *src);
static inline void darray_set_range(struct darray *a,
				    ssize_t i, ssize_t n, const void *src);

/* informative */
#define darray_front(a)          (darray_at(a, 0))
#define darray_set_front(a, val) (darray_set(a, 0, val))
#define darray_back(a)           (darray_at(a, darray_size(a) - 1))
#define darray_set_back(a, val)  (darray_set(a, darray_size(a) - 1, val))

static inline ssize_t darray_size(const struct darray *a);
static inline bool darray_empty(const struct darray *a);
static inline ssize_t darray_capacity(const struct darray *a);
static inline size_t darray_elt_size(const struct darray *a);
static inline ssize_t darray_max_size(const struct darray *a);

/* standard operations */

// insert operations return pionters to the newly created space
void *darray_insert(struct darray *a, ssize_t i, const void *val);
void *darray_insert_all(struct darray *a,
			ssize_t i, const void *ptr, ssize_t n);
void *darray_insert_repeat(struct darray *a,
			   ssize_t i, ssize_t n, const void *val);

void *darray_push_back(struct darray *a, const void *val);
void darray_erase(struct darray *a, ssize_t i);
void darray_erase_range(struct darray *a, ssize_t i, ssize_t n);
void darray_pop_back(struct darray *a);
void darray_clear(struct darray *a);

/* allocation/size modification */
bool darray_reserve(struct darray *a, ssize_t n);
bool darray_resize(struct darray *a, ssize_t n);
bool darray_resize_with(struct darray *a, ssize_t n, const void *val);

/* searching, sorting */
bool darray_contains(const struct darray *a, const void *key, equals_fn equal);
void *darray_find(const struct darray *a, const void *key, equals_fn equal);
ssize_t darray_find_index(const struct darray *a,
			  const void *key, equals_fn equal);
void *darray_find_last(const struct darray *a,
		       const void *key, equals_fn equal);
ssize_t darray_find_last_index(const struct darray *a,
			       const void *key, equals_fn equal);
ssize_t darray_binary_search(const struct darray *a,
			     const void *key, compare_fn compar);
void darray_sort(struct darray *a, compare_fn compar);

/* swap, reverse */
void darray_swap(struct darray *a, ssize_t i, ssize_t j);
void darray_reverse(struct darray *a);

/* inline function definitions */
ssize_t darray_size(const struct darray *a)
{
	return a->size;
}

bool darray_empty(const struct darray *a)
{
	return a->size == 0;
}

ssize_t darray_capacity(const struct darray *a)
{
	return array_size(&a->array);
}

size_t darray_elt_size(const struct darray *a)
{
	return array_elt_size(&a->array);
}

ssize_t darray_max_size(const struct darray *a)
{
	return array_max_size(&a->array);
}

void *darray_at(const struct darray *a, ssize_t i)
{
	assert(0 <= i && i < darray_size(a));
	return array_at(&a->array, i);
}

void darray_set(struct darray *a, ssize_t i, const void *src)
{
	assert(0 <= i && i < darray_size(a));
	array_set(&a->array, i, src);
}

void darray_set_range(struct darray *a, ssize_t i, ssize_t n, const void *src)
{
	assert(0 <= i && i <= darray_size(a) - n);
	array_set_range(&a->array, i, n, src);
}

#endif /* _DARRAY_H */
