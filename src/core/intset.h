#ifndef _INTSET_H
#define _INTSET_H

#include <stddef.h>
#include "array.h"

struct intset {
	struct array keys;
};

struct intset_pos {
	intptr_t key;
	ssize_t index;
};

struct intset_iter {
	const struct intset *set;
	ssize_t index;
};

#define INTSET_KEY(it) (*(intptr_t *)array_item(&(it).set->keys, (it).index))
#define INTSET_FOREACH(it, set) \
	for ((it) = intset_iter_make(set); intset_iter_advance(&(it));)

/* create, destroy */
void intset_init(struct intset *s);
void intset_init_copy(struct intset *s, const struct intset *src);
void intset_assign_copy(struct intset *s, const struct intset *src);
void intset_deinit(struct intset *s);

/* properties */
static inline ssize_t intset_count(const struct intset *s);
static inline intptr_t *intset_items(const struct intset *s);

/* methods */
bool intset_add(struct intset *s, intptr_t key);
bool intset_contains(const struct intset *s, intptr_t key);
void intset_copy_to(const struct intset *s, intptr_t *dst);
void intset_clear(struct intset *s);
bool intset_remove(struct intset *s, intptr_t key);

/* position-based operations */
bool intset_find(const struct intset *s, intptr_t key, struct intset_pos *pos);
void intset_insert(struct intset *s, const struct intset_pos *pos);
void intset_remove_at(struct intset *s, const struct intset_pos *pos);

/* iteration */
struct intset_iter intset_iter_make(const struct intset *s);
void intset_iter_reset(struct intset_iter *it);
bool intset_iter_advance(struct intset_iter *it);

/* inline function definitions */
ssize_t intset_count(const struct intset *s)
{
	assert(s);
	return array_count(&s->keys);
}

static inline intptr_t *intset_items(const struct intset *s)
{
	assert(s);
	return array_to_ptr(&s->keys);
}

#endif /* _INTSET_H */
