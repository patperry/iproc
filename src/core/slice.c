#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "timsort.h"
#include "util.h"
#include "slice.h"

void slice_reverse(struct slice *a)
{
	assert(a);
	if (!slice_count(a))
		return;

	memory_reverse(slice_item(a, 0), slice_count(a), slice_elt_size(a));
}

void slice_copy_to(const struct slice *a, void *dst)
{
	assert(a);
	assert(dst || slice_count(a) == 0);

	if (!slice_count(a))
		return;

	memory_copy_to(slice_item(a, 0), slice_count(a), dst,
		       slice_elt_size(a));
}

void slice_fill(struct slice *a, const void *val)
{
	assert(a);

	if (!slice_count(a))
		return;

	memory_fill(slice_item(a, 0), slice_count(a), val, slice_elt_size(a));
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
			    slice_elt_size(a));
}

ssize_t slice_find_index(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!slice_count(a))
		return -1;

	return forward_find_index(slice_item(a, 0), slice_count(a), match,
				  udata, slice_elt_size(a));
}

void *slice_find_last(const struct slice *a, predicate_fn match, void *udata)
{
	assert(a);
	assert(match);

	if (!slice_count(a))
		return NULL;

	return reverse_find(slice_item(a, 0), slice_count(a), match, udata,
			    slice_elt_size(a));
}

ssize_t slice_find_last_index(const struct slice *a, predicate_fn match,
			      void *udata)
{
	assert(a);
	assert(match);

	if (!slice_count(a))
		return -1;

	return reverse_find_index(slice_item(a, 0), slice_count(a), match,
				  udata, slice_elt_size(a));
}

ssize_t slice_binary_search(const struct slice *a, const void *key,
			    compare_fn compar)
{
	assert(a);
	assert(compar);

	if (!slice_count(a))
		return ~((ssize_t)0);

	return binary_search(slice_item(a, 0), slice_count(a), key, compar,
			     slice_elt_size(a));
}

void slice_sort(struct slice *a, compare_fn compar)
{
	assert(a);
	assert(compar);

	if (!slice_count(a))
		return;

	int err = timsort(slice_item(a, 0), slice_count(a), slice_elt_size(a),
			  compar);
	assert(!err);		// TODO: better error handling
}
