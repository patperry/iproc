#ifndef _LIST_H
#define _LIST_H

/* Dynamic Array type
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include <string.h>		// memcpy
#include "array.h"

struct list {
	struct array array;
	ssize_t count;
};

/* create, destroy */
struct list *list_init(struct list *l, size_t elt_size);
bool list_init_copy(struct list *l, const struct list *src);
void list_deinit(struct list *l);

/* assign */
bool list_assign_array(struct list *l, const void *ptr, ssize_t n);
bool list_assign_copy(struct list *l, const struct list *src);

/* properties */
static inline ssize_t list_capacity(const struct list *a);
bool list_set_capacity(struct list *a, ssize_t n);
static inline ssize_t list_count(const struct list *a);
static inline void *list_item(const struct list *a, ssize_t i);
static inline void list_set_item(struct list *l, ssize_t i, const void *val);
static inline size_t list_elt_size(const struct list *a);



void list_copy_to(const struct list *l, void *dst);
void list_get_range(const struct list *l, ssize_t i, ssize_t n,
		    void *dst);


/* standard operations */

// insert operations return pionters to the newly created space
void *list_add(struct list *l, const void *val);
void *list_add_range(struct list *l, const void *ptr, ssize_t n);
void list_clear(struct list *l);
void *list_insert(struct list *l, ssize_t i, const void *val);
void *list_insert_range(struct list *a,
			ssize_t i, const void *vals, ssize_t n);
void list_remove_at(struct list *l, ssize_t i);
void list_remove_range(struct list *l, ssize_t i, ssize_t n);




/* searching, sorting */
bool list_exists(const struct list *l, predicate_fn match, void *udata);
void *list_find(const struct list *l, predicate_fn match, void *udata);
ssize_t list_find_index(const struct list *l, predicate_fn match, void *udata);
void *list_find_last(const struct list *l, predicate_fn match, void *udata);
ssize_t list_find_last_index(const struct list *l, predicate_fn match, void *udata);
ssize_t list_binary_search(const struct list *l,
			     const void *key, compare_fn compar);
void list_sort(struct list *l, compare_fn compar);

/* reverse */
void list_reverse(struct list *l);

/* inline function definitions */
ssize_t list_count(const struct list *a)
{
	return a->count;
}

ssize_t list_capacity(const struct list *a)
{
	return array_size(&a->array);
}

size_t list_elt_size(const struct list *a)
{
	return array_elt_size(&a->array);
}

void *list_item(const struct list *a, ssize_t i)
{
	assert(0 <= i && i < list_count(a));
	return array_item(&a->array, i);
}

void list_set_item(struct list *l, ssize_t i, const void *val)
{
	assert(0 <= i && i < list_count(l));
	array_set(&l->array, i, val);
}

#endif /* _LIST_H */
