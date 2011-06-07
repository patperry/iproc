#ifndef _UTIL_H
#define _UTIL_H

/* Utility macros and functions.
 *
 * Define ssize_t before including this file.
 */

#include <stddef.h>
#include "compare.h"
#include "delegate.h"

#define MAX(x,y) ((y) > (x) ? (y) : (x))
#define MIN(x,y) ((y) < (x) ? (y) : (x))

#define SWAP(x,y,type) do { type t = (x); (x) = (y); (y) = t; } while (0)

#define container_of(ptr, type, member) \
	((type *)((char *)(ptr) - offsetof(type, member)))

#ifndef alignof
# ifdef HAVE_ALIGNOF
#  define alignof __alignof__
# else
#  define alignof(type) \
          ((sizeof(type) > 1)? offsetof(struct { char c; type x; }, x) : 1)
# endif
#endif

/**
 * ALIGN_PTR:
 * @p: a pointer
 * @a: an alignment (typically a small power of 2)
 *
 * Advance the pointer @p by the minimum number of bytes necessary so that
 * the result is divisible by @a.
 */
#define ALIGN_PTR(p,a) \
	((void *) \
	(((char *)(p) + (size_t)(a) - 1) \
	- \
	((uintptr_t)(((char *)(p) + (size_t)(a) - 1)) \
	& \
	((size_t)(a)-1))))

void *memory_fill(void *begin, ssize_t size, const void *val, size_t elt_size);

void *memory_copy_to(const void *src, ssize_t size, void *dst, size_t elt_size);

void memory_swap(void *val1, void *val2, size_t elt_size);
void *memory_reverse(void *begin, ssize_t size, size_t elt_size);
bool memory_overlaps(const void *ptr1, ssize_t n1,
		     const void *ptr2, ssize_t n2, size_t elt_size);

void *forward_find(const void *begin,
		   ssize_t size,
		   predicate_fn match, void *udata, size_t elt_size);
ssize_t forward_find_index(const void *begin,
			   ssize_t size,
			   predicate_fn match, void *udata, size_t elt_size);
void *reverse_find(const void *begin,
		   ssize_t size,
		   predicate_fn match, void *udata, size_t elt_size);
ssize_t reverse_find_index(const void *begin,
			   ssize_t size,
			   predicate_fn match, void *udata, size_t elt_size);
void *sorted_find(const void *begin,
		  ssize_t size,
		  const void *key, compare_fn compar, size_t elt_size);
ssize_t binary_search(const void *begin, ssize_t size, const void *key,
		      compare_fn compar, size_t elt_size);

void *xcalloc(size_t count, size_t size);
void xfree(void *ptr);
void *xmalloc(size_t size);
void *xrealloc(void *ptr, size_t size);

typedef void (*destroy_fn) (void *val);

#endif /* _UTIL_H */
