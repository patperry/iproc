#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "timsort.h"
#include "util.h"
#include "array.h"

#define INITIAL_CAPACITY   1
#define MIN_CAPACITY_DELTA 4

// 0, 1, 5, 11, 20, 34, 55, 86, 133, 203, 308, ...
static void array_grow(struct array *a, ssize_t delta)
{
	assert(a);
	assert(delta >= 0);
	assert(delta <= SSIZE_MAX - array_count(a));

	ssize_t nmax = SSIZE_MAX;
	ssize_t nmin = array_count(a) + delta;
	ssize_t n = array_capacity(a);
	ssize_t inc;

	if (n >= nmin)
		return;
	
	while (n < nmin && n < nmax) {
		inc = n ? (n >> 1) + MIN_CAPACITY_DELTA	// grow by roughly 1.5
		    : INITIAL_CAPACITY;
		n = (n <= nmax - inc) ? n + inc : nmax;
	}
	array_set_capacity(a, n);
}

void array_init(struct array *a, size_t elt_size)
{
	assert(a);

	a->data = NULL;
	a->count = 0;
	a->capacity = 0;
	a->elt_size = elt_size;
}

void array_init_copy(struct array *a, const struct array *src)
{
	assert(a);
	assert(src);

	array_init(a, array_elt_size(src));
	array_assign_copy(a, src);
}

void array_assign_copy(struct array *a, const struct array *src)
{
	assert(a);
	assert(src);

	ssize_t n = array_count(src);

	array_clear(a);
	array_set_capacity(a, n);
	if (n > 0)
		array_add_range(a, array_item(src, 0), n);
}

void array_deinit(struct array *a)
{
	assert(a);
	xfree(a->data);
}

void array_set_capacity(struct array *a, ssize_t n)
{
	assert(a);
	assert(n >= array_count(a));

	void *data1 = xrealloc(a->data, n * array_elt_size(a));
	a->data = data1;
	a->capacity = n;
}

void *array_add(struct array *a, const void *val)
{
	assert(a);

	return array_insert(a, array_count(a), val);
}

void *array_add_range(struct array *a, const void *vals, ssize_t n)
{
	assert(a);
	assert(n >= 0);
	assert(n <= SSIZE_MAX - array_count(a));

	return array_insert_range(a, array_count(a), vals, n);
}

ssize_t array_binary_search(const struct array *a, const void *key,
			    compare_fn compar)
{
	assert(a);
	assert(compar);

	if (!array_count(a))
		return ~((ssize_t)0);

	return binary_search(array_item(a, 0), array_count(a), key, compar,
			     array_elt_size(a));
}

void array_clear(struct array *a)
{
	assert(a);
	a->count = 0;
}

/* UNNECESSARY contains */
/* UNNECESSARY convert_all */

void array_copy_to(const struct array *a, void *dst)
{
	assert(a);
	assert(dst || array_count(a) == 0);

	if (!array_count(a))
		return;

	memory_copy_to(array_item(a, 0), array_count(a), dst,
		       array_elt_size(a));
}

void array_copy_range_to(const struct array *a, ssize_t i, ssize_t n, void *dst)
{
	assert(0 <= i && i <= array_count(a) - n);
	if (n == 0)
		return;

	memcpy(dst, array_item(a, i), n * array_elt_size(a));
}

/* MISSING equals */

bool array_exists(const struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	return array_find(a, match, udata);
}

void *array_find(const struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!array_count(a))
		return NULL;

	return forward_find(array_item(a, 0), array_count(a), match, udata,
			    array_elt_size(a));

}

ssize_t array_find_all(const struct array *a, predicate_fn match, void *udata,
		       void *dst, ssize_t n)
{
	assert(a);
	assert(match);
	assert(dst || n == 0);

	size_t elt_size = array_elt_size(a);
	ssize_t c = array_count(a);
	ssize_t m = 0;
	ssize_t i;
	const void *val;

	for (i = 0; i < c; i++) {
		val = array_item(a, i);
		if (match(val, udata)) {
			m++;
			if (n > 0) {
				memcpy(dst, val, elt_size);
				dst = (char *)dst + elt_size;
				n--;
			}
		}
	}

	return m;
}

ssize_t array_find_index(const struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!array_count(a))
		return -1;

	return forward_find_index(array_item(a, 0), array_count(a), match,
				  udata, array_elt_size(a));
}

void *array_find_last(const struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!array_count(a))
		return NULL;

	return reverse_find(array_item(a, 0), array_count(a), match, udata,
			    array_elt_size(a));
}

ssize_t array_find_last_index(const struct array *a, predicate_fn match,
			      void *udata)
{
	assert(a);
	assert(match);

	if (!array_count(a))
		return -1;

	return reverse_find_index(array_item(a, 0), array_count(a), match,
				  udata, array_elt_size(a));
}

/* UNNECESSARY for_each */
/* UNNECESSARY get_enumerator */
/* MISSING get_hash_code */
/* UNNECESSARY get_range */
/* UNNECESSARY index_of */

void *array_insert(struct array *a, ssize_t i, const void *val)
{
	assert(a);
	assert(i >= 0);
	assert(i <= array_count(a));

	return array_insert_range(a, i, val, 1);
}

void *array_insert_range(struct array *l, ssize_t i, const void *vals,
			 ssize_t n)
{
	assert(l);
	assert(i >= 0);
	assert(i <= array_count(l));
	assert(n >= 0);
	assert(n <= SSIZE_MAX - array_count(l));

	char *res;
	size_t elt_size = array_elt_size(l);
	ssize_t c0 = array_count(l);

	if (n == 0)
		return NULL;

	array_grow(l, n);
	l->count += n;

	res = array_item(l, i);
	memmove(res + n * elt_size, res, (c0 - i) * elt_size);

	if (vals) {
		memcpy(res, vals, n * elt_size);
	} else {
		memset(res, 0, n * elt_size);
	}

	return res;
}

/* UNNECESSARY last_index_of */

bool array_remove(struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	ssize_t i = array_find_index(a, match, udata);

	if (i >= 0) {
		array_remove_at(a, i);
		return true;
	} else {
		return false;
	}
}

ssize_t array_remove_all(struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	ssize_t r = 0;
	ssize_t n = array_count(a);
	ssize_t i;
	const void *val;

	for (i = n; i > 0; i--) {
		val = array_item(a, i - 1);
		if (match(val, udata)) {
			array_remove_at(a, i - 1);
			r++;
		}
	}

	return r;
}

void array_remove_at(struct array *a, ssize_t i)
{
	assert(a);
	assert(i >= 0);
	assert(i < array_count(a));

	array_remove_range(a, i, 1);
}

void array_remove_range(struct array *a, ssize_t i, ssize_t n)
{
	assert(a);
	assert(0 <= i && i <= array_count(a) - n);
	assert(0 <= n && n <= array_count(a));

	if (n == 0)
		return;

	size_t elt_size = array_elt_size(a);
	ssize_t size0 = array_count(a);
	ssize_t size = size0 - n;
	void *dst = array_item(a, i);
	void *src = (char *)dst + n * elt_size;

	memmove(dst, src, (size0 - (i + n)) * elt_size);
	a->count = size;
}

void array_reverse(struct array *a)
{
	assert(a);

	if (array_count(a))
		return;

	memory_reverse(array_item(a, 0), array_count(a), array_elt_size(a));
}

void array_sort(struct array *a, compare_fn compar)
{
	assert(a);
	assert(compar);

	int err = timsort(array_item(a, 0), array_count(a), array_elt_size(a),
			  compar);
	assert(!err);
}

/* UNNECESSARY to_array */
/* UNNECESSARY to_string */

void array_trim_excess(struct array *a)
{
	assert(a);
	array_set_capacity(a, array_count(a));
}

bool array_true_for_all(const struct array *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	ssize_t n = array_count(a);
	ssize_t i;
	const void *val;

	for (i = 0; i < n; i++) {
		val = array_item(a, i);
		if (!match(val, udata))
			return false;
	}

	return true;
}
