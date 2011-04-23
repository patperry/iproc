#ifndef _intmap_H
#define _intmap_H

/* 
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include "hashset.h"
#include "util.h"

struct intmap {
	struct hashset pairs;
	size_t elt_size;
	size_t elt_align;
	size_t val_offset;
};

struct intmap_pos {
	intptr_t key;
	struct hashset_pos pairs_pos;
};

struct intmap_iter {
	struct hashset_iter pairs_it;
};

/* create, destroy */
bool intmap_init(struct intmap *m, size_t elt_size, size_t elt_align);
bool intmap_init_copy(struct intmap *m, const struct intmap *src);
void intmap_deinit(struct intmap *m);

/* assign, copy, clear */
bool intmap_assign_copy(struct intmap *m, const struct intmap *src);
void *intmap_copy_vals_to(const struct intmap *m, void *dst);
intptr_t *intmap_copy_keys_to(const struct intmap *m, intptr_t *dst);
void intmap_clear(struct intmap *m);

/* informative */
bool intmap_empty(const struct intmap *m);
ssize_t intmap_size(const struct intmap *m);
ssize_t intmap_max_size(const struct intmap *m);
size_t intmap_elt_size(const struct intmap *m);
size_t intmap_elt_align(const struct intmap *m);

bool intmap_contains(const struct intmap *m, intptr_t key);
void *intmap_lookup(const struct intmap *m, intptr_t key);
void *intmap_lookup_with(const struct intmap *m, intptr_t key,
			 const void *val0);

/* modification */
bool intmap_add(struct intmap *m, intptr_t key, const void *val);
ssize_t intmap_add_all(struct intmap *m, const intptr_t *keys,
		       const void *vals, ssize_t n);
void intmap_remove(struct intmap *m, intptr_t key);
void intmap_remove_all(struct intmap *m, const intptr_t *keys, ssize_t n);

/* position-based operations */
void *intmap_find(const struct intmap *m, intptr_t key, struct intmap_pos *pos);
bool intmap_insert(struct intmap *m, struct intmap_pos *pos, const void *val);
void intmap_replace(struct intmap *m, struct intmap_pos *pos, const void *val);
void intmap_erase(struct intmap *m, struct intmap_pos *pos);

/* iteration */
void intmap_iter_init(const struct intmap *m, struct intmap_iter *it);
void intmap_iter_deinit(const struct intmap *m, struct intmap_iter *it);
void intmap_iter_reset(const struct intmap *m, struct intmap_iter *it);
bool intmap_iter_advance(const struct intmap *m, struct intmap_iter *it);
void *intmap_iter_current(const struct intmap *m,
			  const struct intmap_iter *it);
intptr_t intmap_iter_current_key(const struct intmap *m,
				 const struct intmap_iter *it);

#endif /* _intmap_H */
