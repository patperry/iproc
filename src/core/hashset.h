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
	ssize_t enlarge_threshold;  // (table size) * enlarge_factor
	ssize_t shrink_threshold;   // (table size) * shrink_factor
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

bool hashset_init(struct hashset *s, hash_fn hash, equals_fn equals, size_t elt_size);
bool hashset_init_copy(struct hashset *s, const struct hashset *src);
void hashset_deinit(struct hashset *s);

/* assign, copy, clear */
bool hashset_assign_copy(struct hashset *s, const struct hashset *src);
void hashset_clear(struct hashset *s);

/* informative */
ssize_t hashset_size(const struct hashset *s);
bool hashset_empty(const struct hashset *s);
ssize_t hashset_max_size(const struct hashset *s);
size_t hashset_elt_size(const struct hashset *s);
static inline uint32_t hashset_hash(const struct hashset *s, const void *val);
static inline bool hashset_equals(const struct hashset *s, const void *val1,
				  const void *val2);

/* modification */
bool hashset_contains(const struct hashset *s, const void *key);
const void *hashset_lookup(const struct hashset *s, const void *key);
const void *hashset_lookup_with(const struct hashset *s, const void *key,
				const void *val0);
bool hashset_add(struct hashset *s, const void *val);
ssize_t hashset_add_all(struct hashset *s, const void *vals, ssize_t n);
void hashset_remove(struct hashset *s, const void *key);
void hashset_remove_all(struct hashset *s, const void *keys, ssize_t n);

/* position-based operations */
const void *hashset_find(const struct hashset *s, const void *key,
			 struct hashset_pos *pos);
bool hashset_insert(struct hashset *s, struct hashset_pos *pos,
		    const void *val);
void hashset_replace(struct hashset *s, struct hashset_pos *pos,
		     const void *val);
void hashset_erase(struct hashset *s, struct hashset_pos *pos);

/* iteration */
void hashset_iter_init(const struct hashset *s, struct hashset_iter *it);
void hashset_iter_deinit(const struct hashset *s, struct hashset_iter *it);
void hashset_iter_reset(const struct hashset *s, struct hashset_iter *it);
bool hashset_iter_advance(const struct hashset *s, struct hashset_iter *it);
const void * hashset_iter_current(const struct hashset *s,
				  const struct hashset_iter *it);


static inline uint32_t hashset_hash(const struct hashset *s, const void *val)
{
	// if (!s->hash) { return memory_hash(val, hashset_elt_size(s)); }
	return s->hash(val);
}

static inline bool hashset_equals(const struct hashset *s, const void *val1,
				 const void *val2)
{
	// if (!s->equals) { return (memcmp(val1, val2, hashset_elt_size(s)) == 0); }
	return s->equals(val1, val2);
}

#endif /* _HASHSET_H */
