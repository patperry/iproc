#ifndef _SHASHSET_H
#define _SHASHSET_H

#include <stddef.h>

#include "compare.h"
#include "hash.h"
#include "sparsetable.h"

struct shashset {
	struct sparsetable table;
	hash_fn hash;
	equals_fn equals;
	ssize_t enlarge_threshold;	// (table size) * enlarge_factor
};

struct shashset_pos {
	uint32_t hash;
	struct sparsetable_pos insert;
	struct sparsetable_pos existing;
	bool has_existing;
	bool has_insert;
};

struct shashset_iter {
	struct sparsetable_iter table_it;
};
#define SHASHSET_KEY(it) SPARSETABLE_VAL((it).table_it)
#define SHASHSET_FOREACH(it, set) \
	for ((it) = shashset_iter_make(set); shashset_iter_advance(&(it));)

/* create, destroy */
void shashset_init(struct shashset *s, hash_fn hash, equals_fn equals,
		  size_t elt_size);
void shashset_init_copy(struct shashset *s, const struct shashset *src);
void shashset_assign_copy(struct shashset *s, const struct shashset *src);
void shashset_deinit(struct shashset *s);

/* properties */
static inline ssize_t shashset_count(const struct shashset *s);
static inline size_t shashset_elt_size(const struct shashset *s);
static inline bool shashset_equals(const struct shashset *s, const void *key1,
				  const void *key2);
static inline uint32_t shashset_hash(const struct shashset *s, const void *key);

void *shashset_item(const struct shashset *s, const void *key);
void *shashset_set_item(struct shashset *s, const void *key);

/* methods */
void *shashset_add(struct shashset *s, const void *key);
void shashset_clear(struct shashset *s);
bool shashset_contains(const struct shashset *s, const void *key);
bool shashset_remove(struct shashset *s, const void *key);
void shashset_trim_excess(struct shashset *s);

/* position-based operations */
void *shashset_find(const struct shashset *s, const void *key,
		   struct shashset_pos *pos);
void *shashset_insert(struct shashset *s, struct shashset_pos *pos,
		     const void *key);
void shashset_remove_at(struct shashset *s, struct shashset_pos *pos);

/* iteration */
struct shashset_iter shashset_iter_make(const struct shashset *s);
void shashset_iter_reset(struct shashset_iter *it);
bool shashset_iter_advance(struct shashset_iter *it);

/* static method definitions */
ssize_t shashset_count(const struct shashset *s)
{
	return sparsetable_count(&s->table);
}

size_t shashset_elt_size(const struct shashset *s)
{
	return sparsetable_elt_size(&s->table);
}

static inline bool shashset_equals(const struct shashset *s, const void *key1,
				  const void *key2)
{
	return s->equals(key1, key2);
}

static inline uint32_t shashset_hash(const struct shashset *s, const void *key)
{
	return s->hash(key);
}

#endif /* _SHASHSET_H */
