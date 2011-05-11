#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "util.h"
#include "slice.h"

struct slice slice_make(const void *ptr, ssize_t n, size_t elt_size)
{
	struct slice a;
	slice_init(&a, ptr, n, elt_size);
	return a;
}

void slice_init(struct slice *a, const void *data, ssize_t size, size_t elt_size)
{
	assert(a);
	assert(data || size == 0);
	assert(elt_size >= 0);
	assert(size <= SSIZE_MAX / MAX(1, elt_size));

	a->data = (void *)data;
	a->size = size;
	a->elt_size = elt_size;
}

void slice_reverse(struct slice *a)
{
	assert(a);
	if (!slice_count(a))
		return;

	memory_reverse(slice_item(a, 0), slice_count(a), array_elt_size(a));
}

void slice_assign_copy(struct slice *a, const struct slice *src)
{
	assert(a);
	assert(src);
	assert(array_elt_size(a) == array_elt_size(src));
	assert(slice_count(a) == slice_count(src));

	if (!slice_count(a))
		return;

	slice_copy_to(src, slice_item(a, 0));
}

void slice_copy_to(const struct slice *a, void *dst)
{
	assert(a);
	assert(dst || slice_count(a) == 0);

	if (!slice_count(a))
		return;

	memory_copy_to(slice_item(a, 0), slice_count(a), dst,
		       array_elt_size(a));
}

void slice_fill(struct slice *a, const void *val)
{
	assert(a);

	if (!slice_count(a))
		return;

	memory_fill(slice_item(a, 0), slice_count(a), val, array_elt_size(a));
}

bool slice_exists(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);
	return slice_find(a, match, udata);
}

void *slice_find(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!slice_count(a))
		return NULL;

	return forward_find(slice_item(a, 0), slice_count(a), match, udata,
			    array_elt_size(a));
}

ssize_t slice_find_index(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!slice_count(a))
		return -1;

	return forward_find_index(slice_item(a, 0), slice_count(a), match, udata,
				  array_elt_size(a));
}

void *slice_find_last(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);
	
	if (!slice_count(a))
		return NULL;
	
	return reverse_find(slice_item(a, 0), slice_count(a), match, udata,
			    array_elt_size(a));
}

ssize_t slice_find_last_index(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!slice_count(a))
		return -1;

	return reverse_find_index(slice_item(a, 0), slice_count(a), match, udata,
				  array_elt_size(a));
}

ssize_t slice_binary_search(const struct slice *a, const void *key,
			    compare_fn compar)
{
	assert(a);
	assert(compar);

	if (!slice_count(a))
		return ~((ssize_t)0);

	return binary_search(slice_item(a, 0), slice_count(a), key, compar,
			     array_elt_size(a));
}

void slice_sort(struct slice *a, compare_fn compar)
{
	assert(a);
	assert(compar);

	if (!slice_count(a))
		return;

	qsort(slice_item(a, 0), slice_count(a), array_elt_size(a), compar);
}
