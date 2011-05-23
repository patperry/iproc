#include "port.h"

#include <assert.h>
#include <string.h>
#include "util.h"

bool memory_overlaps(const void *ptr1, ssize_t n1,
		     const void *ptr2, ssize_t n2, size_t elt_size)
{
	assert(n1 >= 0);
	assert(n2 >= 0);

	if (n1 == 0 || n2 == 0)
		return false;

	const void *begin1 = ptr1;
	const void *end1 = (char *)begin1 + n1 * elt_size;
	const void *begin2 = ptr2;
	const void *end2 = (char *)begin2 + n2 * elt_size;

	return ((begin1 <= begin2 && begin2 < end1)
		|| (begin2 <= begin1 && begin1 < end2));
}

void *memory_fill(void *begin, ssize_t size, const void *val, size_t elt_size)
{
	assert(begin || size == 0);
	assert(elt_size > 0);

	size_t nbytes = size * elt_size;
	void *end = (char *)begin + nbytes;
	void *dst;

	if (!val) {
		memset(begin, 0, nbytes);
	} else {
		assert(!(begin <= val && val < end));

		for (dst = begin; dst < end; dst = (char *)dst + elt_size) {
			memcpy(dst, val, elt_size);
		}
	}

	return end;
}

void *memory_copy_to(const void *src, ssize_t size, void *dst, size_t elt_size)
{
	assert(src || size == 0);
	assert(size >= 0);
	assert(dst || size == 0);
	assert(elt_size > 0);

	size_t nbytes = size * elt_size;
	memmove(dst, src, nbytes);	// not memcpy; be careful of aliasing
	return (char *)dst + elt_size;
}

void memory_swap(void *val1, void *val2, size_t elt_size)
{
	assert(val1);
	assert(val2);
	assert(elt_size > 0);
	assert(val1 != val2);

	char tmp[elt_size];

	memcpy(tmp, val1, elt_size);
	memcpy(val1, val2, elt_size);
	memcpy(val2, tmp, elt_size);
}

void *memory_reverse(void *begin, ssize_t size, size_t elt_size)
{
	assert(begin || size == 0);
	assert(size >= 0);
	assert(elt_size > 0);

	ssize_t n = size;
	ssize_t i;

	for (i = 0; i < n / 2; i++) {
		memory_swap((char *)begin + i * elt_size,
			    (char *)begin + (n - i - 1) * elt_size, elt_size);
	}

	return (char *)begin + n * elt_size;
}

void *forward_find(const void *begin,
		   ssize_t size,
		   predicate_fn match, void *udata, size_t elt_size)
{
	assert(size >= 0);
	assert(match);
	assert(elt_size > 0);

	const void *end = (char *)begin + size * elt_size;
	const void *ptr;

	for (ptr = begin; ptr < end; ptr = (char *)ptr + elt_size) {
		if (match(ptr, udata))
			return (void *)ptr;
	}

	return NULL;
}

ssize_t forward_find_index(const void *begin,
			   ssize_t size,
			   predicate_fn match, void *udata, size_t elt_size)
{
	assert(size >= 0);
	assert(match);
	assert(elt_size > 0);

	const void *ptr = forward_find(begin, size, match, udata, elt_size);

	if (ptr) {
		return ((char *)ptr - (char *)begin) / elt_size;
	}

	return -1;
}

void *reverse_find(const void *begin,
		   ssize_t size,
		   predicate_fn match, void *udata, size_t elt_size)
{
	assert(size >= 0);
	assert(match);
	assert(elt_size > 0);

	const void *end = (char *)begin + size * elt_size;
	const void *ptr;

	for (ptr = end; ptr > begin;) {
		ptr = (char *)ptr - elt_size;
		if (match(ptr, udata))
			return (void *)ptr;
	}

	return NULL;
}

ssize_t reverse_find_index(const void *begin,
			   ssize_t size,
			   predicate_fn match, void *udata, size_t elt_size)
{
	assert(size >= 0);
	assert(match);
	assert(elt_size > 0);

	const void *ptr = reverse_find(begin, size, match, udata, elt_size);

	if (ptr) {
		return ((char *)ptr - (char *)begin) / elt_size;
	}

	return -1;
}

void *sorted_find(const void *begin, ssize_t size, const void *key,
		  compare_fn compar, size_t elt_size)
{
	assert(size >= 0);
	assert(compar);
	assert(elt_size > 0);

	ssize_t i = binary_search(begin, size, key, compar, elt_size);

	if (i >= 0) {
		return (char *)begin + i * elt_size;
	}

	return NULL;
}

ssize_t binary_search(const void *begin, ssize_t size, const void *key,
		      compare_fn compar, size_t elt_size)
{
	assert(size >= 0);
	assert(compar);
	assert(elt_size > 0);

	ssize_t left = 0;
	ssize_t right = size;
	ssize_t i;
	const void *ptr;
	int cmp;

	while (left < right) {
		i = left + ((right - left) >> 1);
		ptr = (char *)begin + i * elt_size;
		cmp = compar(ptr, key);

		if (cmp < 0) {	// array[i] < key
			left = i + 1;
		} else if (cmp > 0) {	// array[i] > key
			right = i;
		} else {	// array[i] = key
			return i;
		}
	}

	return ~right;		// left == right, not found
}
