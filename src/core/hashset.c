#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "util.h"
#include "hashset.h"

/* The probing method */
/* #define JUMP_(key, num_probes) (1)          // Linear probing */
#define JUMP_(key, num_probes)    (num_probes)	// Quadratic probing

/* How full we let the table get before we resize, by default.
 * Knuth says .8 is good -- higher causes us to probe too much,
 * though it saves memory.
 */
#define HT_OCCUPANCY_PCT 80	// (out of 100);

/* Minimum size we're willing to let hashtables be.
 * Must be a power of two, and at least 4.
 * Note, however, that for a given hashtable, the initial size is a
 * function of the first constructor arg, and may be >HT_MIN_BUCKETS.
 */
#define HT_MIN_BUCKETS	4

/* By default, if you don't specify a hashtable size at
 * construction-time, we use this size.  Must be a power of two, and
 * at least HT_MIN_BUCKETS.
 */
#define HT_DEFAULT_STARTING_BUCKETS	32

#define HT_ENLARGE_FACTOR	(HT_OCCUPANCY_PCT / 100.0)

/* The number of buckets must be a power of 2.  This is the largest 
 * power of 2 that a ssize_t can hold.
 */
#define HT_MAX_BUCKETS	((ssize_t)(((size_t)SSIZE_MAX + 1) >> 1))
#define HT_MAX_SIZE	((ssize_t)(HT_MAX_BUCKETS * HT_ENLARGE_FACTOR))

/* This is the smallest size a hashtable can be without being too crowded
 * If you like, you can give a min #buckets as well as a min #elts */
static ssize_t min_buckets(ssize_t num_elts, ssize_t min_buckets_wanted)
{
	assert(0 <= num_elts && num_elts <= HT_MAX_SIZE);
	assert(0 <= min_buckets_wanted && min_buckets_wanted <= HT_MAX_BUCKETS);

	double enlarge = HT_ENLARGE_FACTOR;
	ssize_t sz = HT_MIN_BUCKETS;	// min buckets allowed

	while (sz < min_buckets_wanted || num_elts > (ssize_t)(sz * enlarge)) {
		assert(sz * 2 > sz);
		sz *= 2;
	}
	return sz;
}

/* Reset the enlarge threshold */
static void hashset_reset_thresholds(struct hashset *s, ssize_t num_buckets)
{
	s->enlarge_threshold = num_buckets * HT_ENLARGE_FACTOR;
}

static ssize_t hashset_bucket_count(const struct hashset *s)
{
	return sparsetable_size(&s->table);
}

static void hashset_init_sized(struct hashset *s, hash_fn hash,
			       equals_fn equals, ssize_t num_buckets,
			       size_t elt_size)
{
	assert(s);
	assert(hash);
	assert(equals);
	assert(num_buckets >= HT_MIN_BUCKETS);
	assert(elt_size >= 0);

	sparsetable_init(&s->table, num_buckets, elt_size);
	s->hash = hash;
	s->equals = equals;
	hashset_reset_thresholds(s, num_buckets);
}

static void hashset_init_copy_sized(struct hashset *s,
				    const struct hashset *src,
				    ssize_t num_buckets)
{
	assert(s);
	assert(src);
	assert(s != src);
	assert(num_buckets >= HT_MIN_BUCKETS);

	struct hashset_iter it;
	const void *val;

	hashset_init_sized(s, src->hash, src->equals, num_buckets,
			   hashset_elt_size(src));

	HASHSET_FOREACH(it, src) {
		val = HASHSET_VAL(it);
		hashset_add(s, val);
	}
}

static bool hashset_needs_grow_delta(const struct hashset *s, ssize_t delta)
{
	assert(delta >= 0);
	assert(sparsetable_count(&s->table) <= HT_MAX_SIZE - delta);

	if (hashset_bucket_count(s) >= HT_MIN_BUCKETS
	    && sparsetable_count(&s->table) + delta <= s->enlarge_threshold) {
		return false;
	} else {
		return true;
	}
}

/* the implementation is simpler than in sparsehashtable.h since we don't
 * store deleted keys */
static void hashset_grow_delta(struct hashset *s, ssize_t delta)
{
	ssize_t num_nonempty = sparsetable_count(&s->table);
	ssize_t bucket_count = hashset_bucket_count(s);
	ssize_t resize_to = min_buckets(num_nonempty + delta, bucket_count);

	struct hashset copy;
	hashset_init_copy_sized(&copy, s, resize_to);
	hashset_deinit(s);
	*s = copy;
}

void hashset_init(struct hashset *s, hash_fn hash, equals_fn equals,
		  size_t elt_size)
{
	assert(s);
	assert(hash);
	assert(equals);
	assert(elt_size >= 0);

	ssize_t num_buckets = HT_DEFAULT_STARTING_BUCKETS;
	hashset_init_sized(s, hash, equals, num_buckets, elt_size);
}

void hashset_init_copy(struct hashset *s, const struct hashset *src)
{
	assert(s);
	assert(src);
	assert(s != src);

	ssize_t num_buckets = hashset_bucket_count(src);
	hashset_init_copy_sized(s, src, num_buckets);
}

void hashset_assign_copy(struct hashset *s, const struct hashset *src)
{
	assert(s);
	assert(src);

	hashset_deinit(s);
	hashset_init_copy(s, src);
}

void hashset_deinit(struct hashset *s)
{
	assert(s);

	sparsetable_deinit(&s->table);
}

void *hashset_item(const struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct hashset_pos pos;
	return hashset_find(s, key, &pos);
}

void *hashset_set_item(struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct hashset_pos pos;
	void *dst;

	if ((dst = hashset_find(s, key, &pos))) {
		memcpy(dst, key, hashset_elt_size(s));
		return dst;
	} else {
		return hashset_insert(s, &pos, key);
	}
}

void *hashset_add(struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct hashset_pos pos;
	if (hashset_find(s, key, &pos)) {
		return NULL;
	} else {
		return hashset_insert(s, &pos, key);
	}
}

void hashset_clear(struct hashset *s)
{
	assert(s);

	ssize_t num_buckets = HT_DEFAULT_STARTING_BUCKETS;

	sparsetable_clear(&s->table);
	sparsetable_set_size(&s->table, num_buckets);
	hashset_reset_thresholds(s, num_buckets);
}

bool hashset_contains(const struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	return hashset_item(s, key);
}

/* MISSING copy_to */
/* INVALID create_set_comparer */
/* MISSING equals */
/* MISSING except_with */
/* INVALID finalize */
/* UNNECESSARY get_enumerator (iter_make) */
/* MISSING get_hash_code */
/* UNNECESSARY get_object_data */
/* INVALID get_type */
/* MISSING intersect_with */
/* MISSING is_proper_subset_of */
/* MISSING is_proper_superset_of */
/* MISSING is_subset_of */
/* MISSING is_superset_of */
/* UNNESSARY memberwise_clone (assign_copy) */
/* INVALID on_deserialization */
/* MISSING overlaps */

bool hashset_remove(struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	bool found;

	struct hashset_pos pos;
	if ((found = hashset_find(s, key, &pos))) {
		hashset_remove_at(s, &pos);
	}
	assert(!hashset_contains(s, key));
	return found;
}

/* MISSING remove_where */
/* MISSING set_equals */
/* MISSING symmetric_except_with */
/* MISSING to_string */

void hashset_trim_excess(struct hashset *s)
{
	assert(s);

	ssize_t count = sparsetable_count(&s->table);
	ssize_t resize_to = min_buckets(count, 0);

	struct hashset copy;
	hashset_init_copy_sized(&copy, s, resize_to);
	hashset_deinit(s);
	*s = copy;
}

/* MISSING union_with */

void *hashset_find(const struct hashset *s, const void *key,
		   struct hashset_pos *pos)
{
	assert(s);
	assert(key);
	assert(pos);

	const struct sparsetable *table = &s->table;
	const ssize_t bucket_count = sparsetable_size(table);
	ssize_t num_probes = 0;	// how many times we've probed
	const ssize_t bucket_count_minus_one = bucket_count - 1;
	uint32_t hash = hashset_hash(s, key);
	ssize_t bucknum = hash & bucket_count_minus_one;
	struct sparsetable_pos table_pos;
	void *val;
	bool deleted;

	pos->hash = hash;
	pos->has_insert = false;
	pos->has_existing = false;

	for (num_probes = 0; num_probes < bucket_count; num_probes++) {
		val = sparsetable_find(table, bucknum, &table_pos);
		deleted = sparsetable_deleted(table, &table_pos);
		if (!val && !deleted) {	// bucket is empty
			if (!pos->has_insert) {	// found no prior place to insert
				pos->insert = table_pos;
				pos->has_insert = true;
			}
			pos->has_existing = false;
			return NULL;
		} else if (!val && deleted) {	// keep searching, but mark to insert
			if (!pos->has_insert) {
				pos->insert = table_pos;
				pos->has_insert = true;
			}
		} else if (hashset_equals(s, key, val)) {
			pos->existing = table_pos;
			pos->has_existing = true;
			return val;
		}
		bucknum =
		    (bucknum +
		     JUMP_(key, num_probes + 1)) & bucket_count_minus_one;
	}

	return NULL;		// table is full and val is not present
}

void *hashset_insert(struct hashset *s, struct hashset_pos *pos,
		     const void *key)
{
	assert(s);
	assert(pos);
	assert(!pos->has_existing);
	assert(key);
	assert(hashset_hash(s, key) == pos->hash);

	if (hashset_needs_grow_delta(s, 1)) {
		hashset_grow_delta(s, 1);
		hashset_find(s, key, pos);	// need to recompute pos
	}

	assert(pos->has_insert);
	return sparsetable_insert(&s->table, &pos->insert, key);
}

void hashset_remove_at(struct hashset *s, struct hashset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(pos->has_existing);

	sparsetable_remove_at(&s->table, &pos->existing);
}

struct hashset_iter hashset_iter_make(const struct hashset *s)
{
	assert(s);

	struct hashset_iter it;
	it.table_it = sparsetable_iter_make(&s->table);
	return it;
}

void hashset_iter_reset(struct hashset_iter *it)
{
	assert(it);

	sparsetable_iter_reset(&it->table_it);
}

bool hashset_iter_advance(struct hashset_iter *it)
{
	assert(it);

	return sparsetable_iter_advance(&it->table_it);
}
