#ifndef _UTIL_H
#define _UTIL_H

/* Utility macros and functions.
 *
 * Define ssize_t before including this file.
 */

#include <stddef.h>
#include "compare.h"

#define MAX(x,y) ((y) > (x) ? (y) : (x))
#define MIN(x,y) ((y) < (x) ? (y) : (x))

#define SWAP(x,y,type) do { type t = (x); (x) = (y); (y) = t; } while (0)

#define container_of(ptr, type, member) \
	((type *)((char *)(ptr) - offsetof(type, member)))

void *memory_fill(void *begin, ssize_t size, const void *val, size_t elt_size);

void *memory_copy_to(const void *src, ssize_t size, void *dst, size_t elt_size);

void memory_swap(void *val1, void *val2, size_t elt_size);
void *memory_reverse(void *begin, ssize_t size, size_t elt_size);

void *forward_find(const void *begin,
		   ssize_t size,
		   const void *key, equals_fn equal, size_t elt_size);
ssize_t forward_find_index(const void *begin,
			   ssize_t size,
			   const void *key, equals_fn equal, size_t elt_size);
void *reverse_find(const void *begin,
		   ssize_t size,
		   const void *key, equals_fn equal, size_t elt_size);
ssize_t reverse_find_index(const void *begin,
			   ssize_t size,
			   const void *key, equals_fn equal, size_t elt_size);
void *sorted_find(const void *begin,
		  ssize_t size,
		  const void *key, compare_fn compar, size_t elt_size);
ssize_t binary_search(const void *begin,
		      ssize_t size,
		      const void *key, compare_fn compar, size_t elt_size);

typedef void (*destroy_fn) (void *val);

#endif /* _UTIL_H */
