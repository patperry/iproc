#ifndef _HASHSET_H
#define _HASHSET_H

#include <stddef.h>

#include "compare.h"
#include "hash.h"
#include "sparsetable.h"

struct hashset {
	struct sparsetable table;
	hash_fn hash;
	equals_fn equals;
	ssize_t enlarge_threshold;	// (table size) * enlarge_factor
};

struct hashset_pos {
	uint32_t hash;
	struct sparsetable_pos insert;
	struct sparsetable_pos existing;
	bool has_existing;
	bool has_insert;
};

struct hashset_iter {
	struct sparsetable_iter table_it;
};
#define HASHSET_KEY(it) SPARSETABLE_VAL((it).table_it)
#define HASHSET_FOREACH(it, set) \
	for ((it) = hashset_iter_make(set); hashset_iter_advance(&(it));)

/* create, destroy */
void hashset_init(struct hashset *s, hash_fn hash, equals_fn equals,
		  size_t elt_size);
void hashset_init_copy(struct hashset *s, const struct hashset *src);
void hashset_assign_copy(struct hashset *s, const struct hashset *src);
void hashset_deinit(struct hashset *s);

/* properties */
static inline ssize_t hashset_count(const struct hashset *s);
static inline size_t hashset_elt_size(const struct hashset *s);
static inline bool hashset_equals(const struct hashset *s, const void *key1,
				  const void *key2);
static inline uint32_t hashset_hash(const struct hashset *s, const void *key);

void *hashset_item(const struct hashset *s, const void *key);
void *hashset_set_item(struct hashset *s, const void *key);

/* methods */
void *hashset_add(struct hashset *s, const void *key);
void hashset_clear(struct hashset *s);
bool hashset_contains(const struct hashset *s, const void *key);
bool hashset_remove(struct hashset *s, const void *key);
void hashset_trim_excess(struct hashset *s);

/* position-based operations */
void *hashset_find(const struct hashset *s, const void *key,
		   struct hashset_pos *pos);
void *hashset_insert(struct hashset *s, struct hashset_pos *pos,
		     const void *key);
void hashset_remove_at(struct hashset *s, struct hashset_pos *pos);

/* iteration */
struct hashset_iter hashset_iter_make(const struct hashset *s);
void hashset_iter_reset(struct hashset_iter *it);
bool hashset_iter_advance(struct hashset_iter *it);

/* static method definitions */
ssize_t hashset_count(const struct hashset *s)
{
	return sparsetable_count(&s->table);
}

size_t hashset_elt_size(const struct hashset *s)
{
	return sparsetable_elt_size(&s->table);
}

static inline bool hashset_equals(const struct hashset *s, const void *key1,
				  const void *key2)
{
	return s->equals(key1, key2);
}

static inline uint32_t hashset_hash(const struct hashset *s, const void *key)
{
	return s->hash(key);
}

#endif /* _HASHSET_H */
