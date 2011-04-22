#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "util.h"
#include "darray.h"

#define INITIAL_CAPACITY   1
#define MIN_CAPACITY_DELTA 4

// 0, 1, 5, 11, 20, 34, 55, 86, 133, 203, 308, ...
static bool darray_grow(struct darray *a)
{
	assert(a);

	ssize_t nmax = darray_max_size(a);
	ssize_t n0 = darray_capacity(a);
	ssize_t inc = n0 ? (n0 >> 1) + MIN_CAPACITY_DELTA	// grow by roughly 1.5
	    : INITIAL_CAPACITY;
	ssize_t n = (n0 <= nmax - inc) ? n0 + inc : nmax;

	if (n != n0 && darray_reserve(a, n)) {
		return a;
	}

	return NULL;
}

static ssize_t darray_reserve_insert(struct darray *a, ssize_t delta)
{
	assert(a);
	assert(delta >= 0);

	do {
		if (darray_size(a) <= darray_capacity(a) - delta)
			return delta;
	} while (darray_grow(a));

	return 0;
}

struct darray *darray_init(struct darray *a, size_t elt_size)
{
	assert(a);
	assert(elt_size > 0);

	if (array_init(&a->array, 0, elt_size)) {
		a->size = 0;
		return a;
	}

	return NULL;
}

bool darray_init_copy(struct darray * a, const struct darray * src)
{
	assert(a);
	assert(src);

	if (darray_init(a, darray_elt_size(src))) {
		if (darray_assign_copy(a, src)) {
			return a;
		}
		darray_deinit(a);
	}

	return NULL;
}

void darray_deinit(struct darray *a)
{
	assert(a);
	array_deinit(&a->array);
}

bool darray_assign_repeat(struct darray *a, ssize_t n, const void *val)
{
	assert(a);
	assert(n >= 0);
	assert(n <= darray_max_size(a));

	if (darray_reserve(a, n)) {
		darray_clear(a);
		darray_insert_repeat(a, 0, n, val);
		return a;
	}

	return NULL;
}

bool darray_assign(struct darray * a, const void *ptr, ssize_t n)
{
	assert(a);
	assert(n >= 0);
	assert(n <= darray_max_size(a));

	if (darray_reserve(a, n)) {
		darray_clear(a);
		darray_insert_all(a, 0, ptr, n);
		return a;
	}

	return NULL;
}

bool darray_assign_copy(struct darray * a, const struct darray * src)
{
	assert(a);
	assert(src);

	const void *psrc = darray_empty(src) ? NULL : darray_front(src);
	return darray_assign(a, psrc, darray_size(src));
}

void *darray_copy_to(const struct darray *a, void *dst)
{
	assert(a);
	assert(dst || darray_size(a) == 0);

	if (darray_empty(a))
		return dst;
	
	return memory_copy_to(darray_front(a), darray_size(a), dst,
			      darray_elt_size(a));
}

void *darray_copy_range_to(const struct darray *a, ssize_t i, ssize_t n, void *dst)
{
	assert(0 <= i && i <= darray_size(a) - n);
	return array_copy_range_to(&a->array, i, n, dst);
}

void darray_fill(struct darray *a, const void *val)
{
	assert(a);
	darray_fill_range(a, 0, darray_size(a), val);
}

void darray_fill_range(struct darray *a, ssize_t i, ssize_t n, const void *val)
{
	assert(a);
	assert(n >= 0);
	assert(0 <= i && i <= darray_size(a) - n);

	if (n == 0)
		return;
	
	memory_fill(darray_at(a, i), n, val, darray_elt_size(a));
}

static ssize_t darray_insert_space(struct darray *a, ssize_t i, ssize_t n)
{
	ssize_t size0 = a->size;
	ssize_t size = size0 + n;
	size_t elt_size = darray_elt_size(a);
	size_t tail_size = (size0 - i) * elt_size;

	void *src, *dst;
	ssize_t nactual;
	
	if ((nactual = darray_reserve_insert(a, n))) {
		darray_resize_with(a, size, NULL);
		src = darray_at(a, i);	// compute dst after resize in case of realloc
		dst = (char *)src + nactual * elt_size;
		memmove(dst, src, tail_size);
	}

	return nactual;
}

bool darray_insert(struct darray *a, ssize_t i, const void *val)
{
	assert(a);
	assert(darray_size(a) <= darray_max_size(a) - 1);
	assert(i >= 0);
	assert(i <= darray_size(a));

	return darray_insert_repeat(a, i, 1, val);
}

ssize_t darray_insert_repeat(struct darray *a, ssize_t i, ssize_t n,
			   const void *val)
{
	assert(a);
	assert(i >= 0);
	assert(i <= darray_size(a));
	assert(n >= 0);
	assert(n <= darray_max_size(a) - darray_size(a));

	ssize_t nactual;
	
	if ((nactual = darray_insert_space(a, i, n))) {
		if (val)
			array_fill_range(&a->array, i, nactual, val);
	}

	return nactual;
}

ssize_t darray_insert_all(struct darray *a, ssize_t i, const void *ptr, ssize_t n)
{
	assert(a);
	assert(i >= 0);
	assert(i <= darray_size(a));
	assert(ptr || n == 0);
	assert(n >= 0);
	assert(n <= darray_max_size(a) - darray_size(a));

	ssize_t nactual;

	if ((nactual = darray_insert_space(a, i, n)) && ptr) {
		array_set_range(&a->array, i, nactual, ptr);
	}

	return nactual;
}

bool darray_push_back(struct darray *a, const void *val)
{
	assert(a);
	assert(darray_size(a) <= darray_max_size(a) - 1);

	return darray_insert(a, darray_size(a), val);
}

void darray_erase(struct darray *a, ssize_t i)
{
	assert(a);
	assert(i >= 0);
	assert(i < darray_size(a));

	darray_erase_range(a, i, 1);
}

void darray_erase_range(struct darray *a, ssize_t i, ssize_t n)
{
	assert(a);
	assert(i >= 0);
	assert(i <= darray_size(a) - n);
	assert(n >= 0);

	size_t elt_size = darray_elt_size(a);
	ssize_t size = darray_size(a) - n;
	void *dst = darray_at(a, i);
	void *src = (char *)dst + n * elt_size;

	memmove(dst, src, n * elt_size);
	darray_resize(a, size);
}

void darray_pop_back(struct darray *a)
{
	assert(a);
	assert(!darray_empty(a));

	darray_erase(a, darray_size(a) - 1);
}

void darray_clear(struct darray *a)
{
	assert(a);
	a->size = 0;
}

bool darray_reserve(struct darray *a, ssize_t n)
{
	assert(a);
	assert(n <= darray_max_size(a));

	if (darray_capacity(a) >= n
	    || array_reinit(&a->array, n, darray_elt_size(a))) {
		return true;
	}
	return false;
}

bool darray_resize(struct darray * a, ssize_t n)
{
	assert(a);
	assert(n >= 0);
	assert(n <= darray_max_size(a));

	ssize_t n0 = a->size;

	if (darray_resize_with(a, n, NULL)) {
		if (n > n0)
			array_fill_range(&a->array, n0, n - n0, NULL);	// fill with zeros
		return a;
	}

	return NULL;
}

bool darray_resize_with(struct darray * a, ssize_t n, const void *val)
{
	assert(a);
	assert(n >= 0);
	assert(n <= darray_max_size(a));

	ssize_t n0 = a->size;

	if (darray_reserve(a, n)) {
		a->size = n;
		if (n > n0 && val)
			array_fill_range(&a->array, n0, n - n0, val);
		return a;
	}

	return NULL;
}

bool darray_contains(const struct darray * a, const void *key, equals_fn equal)
{
	assert(a);
	assert(equal);

	return darray_find(a, key, equal);
}

void *darray_find(const struct darray *a, const void *key, equals_fn equal)
{
	assert(a);
	assert(equal);

	if (darray_empty(a))
		return NULL;
	
	return forward_find(darray_front(a), darray_size(a), key, equal,
			    darray_elt_size(a));

}

ssize_t darray_find_index(const struct darray * a, const void *key,
			  equals_fn equal)
{
	assert(a);
	assert(equal);

	if (darray_empty(a))
		return -1;
	
	return forward_find_index(darray_front(a), darray_size(a), key, equal,
				  darray_elt_size(a));
}

void *darray_find_last(const struct darray *a, const void *key, equals_fn equal)
{
	assert(a);
	assert(equal);

	if (darray_empty(a))
		return NULL;
	
	return reverse_find(darray_front(a), darray_size(a), key, equal,
			    darray_elt_size(a));
}

ssize_t darray_find_last_index(const struct darray * a, const void *key,
			       equals_fn equal)
{
	assert(a);
	assert(equal);

	if (darray_empty(a))
		return -1;
	
	return reverse_find_index(darray_front(a), darray_size(a), key, equal,
				  darray_elt_size(a));
}

ssize_t darray_binary_search(const struct darray * a, const void *key,
			     compare_fn compar)
{
	assert(a);
	assert(compar);

	if (darray_empty(a))
		return ~((ssize_t)0);
	
	return binary_search(darray_front(a), darray_size(a), key, compar,
			     darray_elt_size(a));
}

void darray_sort(struct darray *a, compare_fn compar)
{
	assert(a);
	assert(compar);

	qsort(darray_front(a), darray_size(a), darray_elt_size(a), compar);
}

void darray_swap(struct darray *a, ssize_t i, ssize_t j)
{
	assert(a);
	assert(0 <= i && i < darray_size(a));
	assert(0 <= j && j < darray_size(a));
	assert(i != j);

	memory_swap(darray_at(a, i), darray_at(a, j), darray_elt_size(a));
}

void darray_reverse(struct darray *a)
{
	assert(a);
	
	if (darray_empty(a))
		return;
	
	memory_reverse(darray_front(a), darray_size(a), darray_elt_size(a));
}
