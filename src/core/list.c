#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "util.h"
#include "list.h"

#define INITIAL_CAPACITY   1
#define MIN_CAPACITY_DELTA 4

// 0, 1, 5, 11, 20, 34, 55, 86, 133, 203, 308, ...
static bool list_grow(struct list *l, ssize_t delta)
{
	assert(l);
	assert(delta >= 0);
	assert(delta <= SSIZE_MAX - list_count(l));

	ssize_t nmax = SSIZE_MAX;
	ssize_t nmin = list_count(l) + delta;
	ssize_t n = list_capacity(l);
	ssize_t inc;
	
	while (n < nmin && n < nmax) {
		inc = n ? (n >> 1) + MIN_CAPACITY_DELTA	// grow by roughly 1.5
			: INITIAL_CAPACITY;
		n = (n <= nmax - inc) ? n + inc : nmax;
	}
	
	return list_set_capacity(l, n);
}

bool list_init(struct list *l, size_t elt_size)
{
	assert(l);
	assert(elt_size >= 0);

	l->data = NULL;
	l->count = 0;
	l->capacity = 0;
	l->elt_size = elt_size;
	return true;
}

bool list_init_copy(struct list *l, const struct list *src)
{
	assert(l);
	assert(src);

	if (list_init(l, list_elt_size(src))) {
		if (list_assign_copy(l, src)) {
			return l;
		}
		list_deinit(l);
	}

	return NULL;
}

void list_deinit(struct list *l)
{
	assert(l);
	free(l->data);
}

bool list_assign_array(struct list *l, const void *ptr, ssize_t n)
{
	assert(l);
	assert(n >= 0);

	if (list_set_capacity(l, n)) {
		list_clear(l);
		if (n > 0)
			list_insert_range(l, 0, ptr, n);
		return l;
	}

	return NULL;
}

bool list_assign_copy(struct list *l, const struct list *src)
{
	assert(l);
	assert(src);

	const void *psrc = list_count(src) ? list_item(src, 0) : NULL;
	return list_assign_array(l, psrc, list_count(src));
}

bool list_set_capacity(struct list *l, ssize_t n)
{
	assert(l);
	assert(n >= list_count(l));
	
	void *data1 = realloc(l->data, n * list_elt_size(l));
	if (data1) {
		l->data = data1;
		l->capacity = n;
		return true;
	}
	
	return false;
}

void *list_add(struct list *l, const void *val)
{
	assert(l);
	
	return list_insert(l, list_count(l), val);
}

void *list_add_range(struct list *l, const void *ptr, ssize_t n)
{
	assert(l);
	assert(n >= 0);
	assert(n <= SSIZE_MAX - list_count(l));
	
	return list_insert_range(l, list_count(l), ptr, n);
}

ssize_t list_binary_search(const struct list *l, const void *key,
			     compare_fn compar)
{
	assert(l);
	assert(compar);
	
	if (!list_count(l))
		return ~((ssize_t)0);
	
	return binary_search(list_item(l, 0), list_count(l), key, compar,
			     list_elt_size(l));
}

void list_clear(struct list *l)
{
	assert(l);
	l->count = 0;
}

/* MISSING contains */
/* MISSING convert_all */

void list_copy_to(const struct list *l, void *dst)
{
	assert(l);
	assert(dst || list_count(l) == 0);
	
	if (!list_count(l))
		return;
	
	memory_copy_to(list_item(l, 0), list_count(l), dst,
		       list_elt_size(l));
}

void list_copy_range_to(const struct list *l, ssize_t i, ssize_t n, void *dst)
{
	assert(0 <= i && i <= list_count(l) - n);
	if (n == 0)
		return;
	
	memcpy(dst, list_item(l, i), n * list_elt_size(l));
}

/* MISSING equals */

bool list_exists(const struct list *l, predicate_fn match, void *udata)
{
	assert(l);
	assert(match);
	
	return list_find(l, match, udata);
}

void *list_find(const struct list *l, predicate_fn match, void *udata)
{
	assert(l);
	assert(match);
	
	if (!list_count(l))
		return NULL;
	
	return forward_find(list_item(l, 0), list_count(l), match, udata,
			    list_elt_size(l));
	
}

/* MISSING find_all */

ssize_t list_find_index(const struct list *l, predicate_fn match, void *udata)
{
	assert(l);
	assert(match);
	
	if (!list_count(l))
		return -1;
	
	return forward_find_index(list_item(l, 0), list_count(l), match, udata,
				  list_elt_size(l));
}

void *list_find_last(const struct list *l, predicate_fn match, void *udata)
{
	assert(l);
	assert(match);
	
	if (!list_count(l))
		return NULL;
	
	return reverse_find(list_item(l, 0), list_count(l), match, udata,
			    list_elt_size(l));
}

ssize_t list_find_last_index(const struct list *l, predicate_fn match, void *udata)
{
	assert(l);
	assert(match);

	if (!list_count(l))
		return -1;
	
	return reverse_find_index(list_item(l, 0), list_count(l), match, udata,
				  list_elt_size(l));
}

/* MISSING for_each */
/* MISSING get_enumerator */
/* MISSING get_hash_code */
/* MISSING get_range */
/* MISSING index_of */

void *list_insert(struct list *l, ssize_t i, const void *val)
{
	assert(l);
	assert(i >= 0);
	assert(i <= list_count(l));

	return list_insert_range(l, i, val, 1);
}

void *list_insert_range(struct list *l, ssize_t i, const void *vals, ssize_t n)
{
	assert(l);
	assert(i >= 0);
	assert(i <= list_count(l));
	assert(n >= 0);
	assert(n <= SSIZE_MAX - list_count(l));

	char *res;
	size_t elt_size = list_elt_size(l);
	ssize_t c0 = list_count(l);
	
	if (n == 0)
		return NULL;
	
	if (!list_grow(l, n))
		return NULL;
	
	l->count += n;
	
	res = list_item(l, i);
	memmove(res + n * elt_size, res, (c0 - i) * elt_size);
	
	if (vals) {
		memcpy(res, vals, n * elt_size);
	} else {
		memset(res, 0, n * elt_size);
	}

	return res;
}

/* MISSING last_index_of */
/* MISSING remove */
/* MISSING remove_all */

void list_remove_at(struct list *l, ssize_t i)
{
	assert(l);
	assert(i >= 0);
	assert(i < list_count(l));
	
	list_remove_range(l, i, 1);
}

void list_remove_range(struct list *l, ssize_t i, ssize_t n)
{
	assert(l);
	assert(i >= 0);
	assert(i <= list_count(l) - n);
	assert(n >= 0);
	
	if (n == 0)
		return;
	
	size_t elt_size = list_elt_size(l);
	ssize_t size = list_count(l) - n;
	void *dst = list_item(l, i);
	void *src = (char *)dst + n * elt_size;
	
	memmove(dst, src, n * elt_size);
	l->count = size;
}

void list_reverse(struct list *l)
{
	assert(l);
	
	if (list_count(l))
		return;
	
	memory_reverse(list_item(l, 0), list_count(l), list_elt_size(l));
}

void list_sort(struct list *l, compare_fn compar)
{
	assert(l);
	assert(compar);
	
	qsort(list_item(l, 0), list_count(l), list_elt_size(l), compar);
}

/* UNNECESSARY to_array */
/* MISSING to_string */
/* MISSING trim_excess */
/* MISSING true_for_all */
