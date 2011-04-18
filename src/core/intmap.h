#ifndef _intmap_H
#define _intmap_H

/* 
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include "darray.h"
#include "util.h"

struct intmap_item {
	intptr_t key;
	intptr_t val;
};

struct intmap {
	struct darray items;
};

struct intmap_pos {
	intptr_t key;
	ssize_t  index;
};

struct intmap_iter {
	ssize_t index;
};


/* create, destroy */
bool intmap_init(struct intmap *m);
bool intmap_init_copy(struct intmap *m, const struct intmap *src);
void intmap_deinit(struct intmap *m);

/* assign, copy, clear */
bool intmap_assign_copy(struct intmap *m, const struct intmap *src);
struct intmap_item *intmap_copy_to(const struct intmap *m,
				   struct intmap_item *dsst);
intptr_t *intmap_copy_vals_to(const struct intmap *m, intptr_t *dst);
intptr_t *intmap_copy_keys_to(const struct intmap *m, intptr_t *dst);
void intmap_clear(struct intmap *m);

/* informative */
bool intmap_empty(const struct intmap *m);
ssize_t intmap_size(const struct intmap *m);
ssize_t intmap_max_size(const struct intmap *m);

bool intmap_contains(const struct intmap *m, intptr_t key);
intptr_t intmap_lookup(const struct intmap *m, intptr_t key);
intptr_t intmap_lookup_with(const struct intmap *m, intptr_t key,
			    intptr_t val0);


/* modification */
bool intmap_add(struct intmap *m, intptr_t key, intptr_t val);
bool intmap_add_all(struct intmap *m, const intptr_t *keys,
		    const intptr_t *vals, ssize_t n);
void intmap_remove(struct intmap *m, intptr_t key);
void intmap_remove_all(struct intmap *m, const intptr_t *keys, ssize_t n);

/* position-based operations */
const intptr_t *intmap_find(const struct intmap *m, intptr_t key,
			    struct intmap_pos *pos);
bool intmap_insert(struct intmap *m, struct intmap_pos *pos,
		   intptr_t val);
void intmap_replace(struct intmap *m, struct intmap_pos *pos, intptr_t val);
void intmap_erase(struct intmap *m, struct intmap_pos *pos);

/* iteration */
void intmap_iter_init(const struct intmap *m, struct intmap_iter *it);
void intmap_iter_deinit(const struct intmap *m, struct intmap_iter *it);
void intmap_iter_reset(const struct intmap *m, struct intmap_iter *it);
bool intmap_iter_advance(const struct intmap *m, struct intmap_iter *it);
struct intmap_item intmap_iter_current(const struct intmap *m,
				       const struct intmap_iter *it);







#endif /* _intmap_H */
