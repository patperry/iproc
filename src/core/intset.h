#ifndef _INTSET_H
#define _INTSET_H

#include <stddef.h>
#include "darray.h"

struct intset {
	struct darray values;
};

struct intset_pos {
	intptr_t value;
	ssize_t  index;
};

struct intset_iter {
	ssize_t index;
};

/* create, destroy */
bool intset_init(struct intset *s);
bool intset_init_copy(struct intset *s, const struct intset *src);
void intset_deinit(struct intset *s);

/* assign, copy, clear */
bool intset_assign(struct intset *s, const intptr_t *vals, ssize_t n);
bool intset_assign_sorted(struct intset *s, const intptr_t *vals, ssize_t n);
bool intset_assign_copy(struct intset *s, const struct intset *src);
intptr_t *intset_copy_to(const struct intset *s, intptr_t *dst);
void intset_clear(struct intset *s);

/* informative */
bool intset_empty(const struct intset *s);
ssize_t intset_size(const struct intset *s);
ssize_t intset_max_size(const struct intset *s);
intptr_t intset_min(const struct intset *s);
intptr_t intset_max(const struct intset *s);
bool intset_contains(const struct intset *s, intptr_t val);

/* modification */
bool intset_add(struct intset *s, intptr_t val);
bool intset_add_all(struct intset *s, const intptr_t *vals, ssize_t n);
void intset_remove(struct intset *s, intptr_t val);
void intset_remove_all(struct intset *s, const intptr_t *vals, ssize_t n);

/* position-based operations */
bool intset_find(const struct intset *s, intptr_t val, struct intset_pos *pos);
bool intset_insert(struct intset *s, const struct intset_pos *pos);
void intset_erase(struct intset *s, const struct intset_pos *pos);

/* iteration */
void intset_iter_init(const struct intset *s, struct intset_iter *it);
void intset_iter_deinit(const struct intset *s, struct intset_iter *it);
void intset_iter_reset(const struct intset *s, struct intset_iter *it);
bool intset_iter_advance(const struct intset *s, struct intset_iter *it);
intptr_t intset_iter_current(const struct intset *s, const struct intset_iter *it);



#endif /* _INTSET_H */
