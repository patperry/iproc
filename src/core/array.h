#ifndef _ARRAY_H
#define _ARRAY_H

/* Dynamic Array type
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include <string.h>		// memcpy
#include "compare.h"
#include "delegate.h"

struct array {
	void *data;
	ssize_t count;
	ssize_t capacity;
	size_t elt_size;
};

/* create, destroy */
void array_init(struct array *a, size_t elt_size);
void array_init_copy(struct array *a, const struct array *src);
void array_assign_copy(struct array *a, const struct array *src);
void array_deinit(struct array *a);

/* properties */
static inline ssize_t array_capacity(const struct array *a);
void array_set_capacity(struct array *a, ssize_t n);
static inline ssize_t array_count(const struct array *a);
static inline void *array_item(const struct array *a, ssize_t i);
static inline void array_set_item(struct array *a, ssize_t i, const void *val);
static inline size_t array_elt_size(const struct array *a);
static inline void *array_to_ptr(const struct array *a);

void array_copy_to(const struct array *a, void *dst);
void array_copy_range_to(const struct array *a, ssize_t i, ssize_t n,
			 void *dst);

/* standard operations */

// insert operations return pointers to the newly created space
void *array_add(struct array *a, const void *val);
void *array_add_range(struct array *a, const void *vals, ssize_t n);
void array_clear(struct array *a);
void *array_insert(struct array *a, ssize_t i, const void *val);
void *array_insert_range(struct array *a,
			 ssize_t i, const void *vals, ssize_t n);
bool array_remove(struct array *a, predicate_fn match, void *udata);
ssize_t array_remove_all(struct array *a, predicate_fn match, void *udata);
void array_remove_at(struct array *a, ssize_t i);
void array_remove_range(struct array *a, ssize_t i, ssize_t n);

/* searching, sorting */
bool array_exists(const struct array *a, predicate_fn match, void *udata);
void *array_find(const struct array *a, predicate_fn match, void *udata);
ssize_t array_find_all(const struct array *a, predicate_fn match, void *udata,
		       void *dst, ssize_t n);
ssize_t array_find_index(const struct array *a, predicate_fn match,
			 void *udata);
void *array_find_last(const struct array *a, predicate_fn match, void *udata);
ssize_t array_find_last_index(const struct array *a, predicate_fn match,
			      void *udata);
bool array_true_for_all(const struct array *a, predicate_fn match, void *udata);

ssize_t array_binary_search(const struct array *a, const void *key,
			    compare_fn compar);
void array_reverse(struct array *a);
void array_sort(struct array *a, compare_fn compar);
void array_trim_excess(struct array *a);

#define ARRAY_FOREACH(val, a) \
	for ((val) = (a)->data; \
	     (val) < (void *)((char *)(a)->data + (a)->count * (a)->elt_size); \
	     (val)++)

/* inline function definitions */
ssize_t array_count(const struct array *a)
{
	return a->count;
}

ssize_t array_capacity(const struct array *a)
{
	return a->capacity;
}

size_t array_elt_size(const struct array *a)
{
	return a->elt_size;
}

void *array_to_ptr(const struct array *a)
{
	return (void *)a->data;
}

void *array_item(const struct array *a, ssize_t i)
{
	assert(0 <= i && i < array_count(a));
	return (char *)a->data + i * array_elt_size(a);
}

void array_set_item(struct array *a, ssize_t i, const void *val)
{
	assert(0 <= i && i < array_count(a));

	void *ptr = array_item(a, i);
	memcpy(ptr, val, array_elt_size(a));
}

#endif /* _ARRAY_H */
