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

/* How empty we let the table get before we resize lower, by default.
 * (0.0 means never resize lower.)
 * It should be less than OCCUPANCY_PCT / 2 or we thrash resizing
 */
#define HT_EMPTY_PCT	(0.4 * HT_OCCUPANCY_PCT)

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
#define HT_SHRINK_FACTOR	(HT_EMPTY_PCT / 100.0)

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

/* Reset the enlarge and shrink thresholds */
static void hashset_reset_thresholds(struct hashset *s, ssize_t num_buckets)
{
	s->enlarge_threshold = num_buckets * HT_ENLARGE_FACTOR;
	s->shrink_threshold = num_buckets * HT_SHRINK_FACTOR;
}

static ssize_t hashset_bucket_count(const struct hashset *s)
{
	return sparsetable_size(&s->table);
}

static bool hashset_init_sized(struct hashset *s, hash_fn hash,
			       equals_fn equals, ssize_t num_buckets,
			       size_t elt_size)
{
	assert(s);
	assert(hash);
	assert(equals);
	assert(num_buckets >= HT_MIN_BUCKETS);
	assert(elt_size >= 0);

	if (sparsetable_init(&s->table, num_buckets, elt_size)) {
		s->hash = hash;
		s->equals = equals;
		hashset_reset_thresholds(s, num_buckets);
		return true;
	}

	return false;
}

bool hashset_init(struct hashset *s, hash_fn hash, equals_fn equals,
		  size_t elt_size)
{
	assert(s);
	assert(hash);
	assert(equals);
	assert(elt_size >= 0);

	ssize_t num_buckets = HT_DEFAULT_STARTING_BUCKETS;
	return hashset_init_sized(s, hash, equals, num_buckets, elt_size);
}

static bool hashset_init_copy_sized(struct hashset *s,
				    const struct hashset *src,
				    ssize_t num_buckets)
{
	assert(s);
	assert(src);
	assert(s != src);
	assert(num_buckets >= HT_MIN_BUCKETS);

	struct hashset_iter it;
	const void *val;
	bool result = true;

	if (!hashset_init_sized(s, src->hash, src->equals, num_buckets,
				hashset_elt_size(src)))
		return false;

	hashset_iter_init(src, &it);
	while (hashset_iter_advance(src, &it)) {
		val = hashset_iter_current(src, &it);
		if (!hashset_add(s, val)) {
			result = false;
			hashset_deinit(s);
			break;
		}

	}
	hashset_iter_deinit(src, &it);
	return result;
}

bool hashset_init_copy(struct hashset *s, const struct hashset *src)
{
	assert(s);
	assert(src);
	assert(s != src);

	ssize_t num_buckets = hashset_bucket_count(src);
	return hashset_init_copy_sized(s, src, num_buckets);
}

void hashset_deinit(struct hashset *s)
{
	assert(s);
	sparsetable_deinit(&s->table);
}

bool hashset_assign_copy(struct hashset *s, const struct hashset *src)
{
	struct hashset tmp;

	if (hashset_init_copy(&tmp, src)) {
		hashset_deinit(s);
		*s = tmp;
		return true;
	}
	return false;
}

void hashset_clear(struct hashset *s)
{
	assert(s);

	ssize_t num_buckets = HT_DEFAULT_STARTING_BUCKETS;

	sparsetable_resize(&s->table, num_buckets);
	sparsetable_clear(&s->table);
	hashset_reset_thresholds(s, num_buckets);
}

ssize_t hashset_size(const struct hashset *s)
{
	return sparsetable_count(&s->table);
}

bool hashset_empty(const struct hashset *s)
{
	return hashset_size(s) == 0;
}

ssize_t hashset_max_size(const struct hashset *s)
{
	return MIN(sparsetable_max_size(&s->table), HT_MAX_SIZE);
}

size_t hashset_elt_size(const struct hashset *s)
{
	return sparsetable_elt_size(&s->table);
}

static bool hashset_needs_shrink(const struct hashset *s)
{
	ssize_t num_remain = sparsetable_count(&s->table);
	ssize_t shrink_threshold = s->shrink_threshold;
	ssize_t bucket_count = hashset_bucket_count(s);

	/* If you construct a hashtable with < HT_DEFAULT_STARTING_BUCKETS,
	 * we'll never shrink until you get relatively big, and we'll never
	 * shrink below HT_DEFAULT_STARTING_BUCKETS.  Otherwise, something
	 * like "dense_hash_set<int> x; x.insert(4); x.erase(4);" will
	 * shrink us down to HT_MIN_BUCKETS buckets, which is too small.
	 */
	return (num_remain < shrink_threshold &&
		bucket_count > HT_DEFAULT_STARTING_BUCKETS);
}

static bool hashset_shrink(struct hashset *s)
{
	assert(hashset_needs_shrink(s));

	ssize_t num_remain = sparsetable_count(&s->table);
	ssize_t bucket_count = hashset_bucket_count(s);
	double shrink_factor = HT_SHRINK_FACTOR;
	ssize_t sz = bucket_count / 2;	// find how much we should shrink

	while (sz > HT_DEFAULT_STARTING_BUCKETS &&
	       num_remain < (ssize_t)(sz * shrink_factor)) {
		sz /= 2;	// stay a power of 2
	}

	struct hashset copy;
	if (hashset_init_copy_sized(&copy, s, sz)) {
		hashset_deinit(s);
		*s = copy;
		return true;
	}

	return false;
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
static bool hashset_grow_delta(struct hashset *s, ssize_t delta)
{
	ssize_t num_nonempty = sparsetable_count(&s->table);
	ssize_t bucket_count = hashset_bucket_count(s);
	ssize_t resize_to = min_buckets(num_nonempty + delta, bucket_count);

	struct hashset copy;
	if (hashset_init_copy_sized(&copy, s, resize_to)) {
		hashset_deinit(s);
		*s = copy;
		return true;
	}

	return false;
}

bool hashset_contains(const struct hashset *s, const void *key)
{
	return hashset_lookup(s, key);
}

const void *hashset_lookup(const struct hashset *s, const void *key)
{
	return hashset_lookup_with(s, key, NULL);
}

const void *hashset_lookup_with(const struct hashset *s, const void *key,
				const void *val0)
{
	struct hashset_pos pos;
	const void *val = hashset_find(s, key, &pos);
	return val ? val : val0;
}

bool hashset_add(struct hashset *s, const void *val)
{
	struct hashset_pos pos;
	if (hashset_find(s, val, &pos)) {
		hashset_replace(s, &pos, val);
		assert(hashset_contains(s, val));
		return true;
	} else {
		return hashset_insert(s, &pos, val);
	}
}

ssize_t hashset_add_all(struct hashset *s, const void *vals, ssize_t n)
{
	size_t elt_size = hashset_elt_size(s);
	ssize_t i;

	for (i = 0; i < n; i++, vals += elt_size) {
		if (!hashset_add(s, vals))
			break;
	}

	return i;
}

void hashset_remove(struct hashset *s, const void *key)
{
	struct hashset_pos pos;
	if (hashset_find(s, key, &pos)) {
		hashset_erase(s, &pos);
	}
	assert(!hashset_contains(s, key));
}

void hashset_remove_all(struct hashset *s, const void *keys, ssize_t n)
{
	size_t elt_size = hashset_elt_size(s);
	ssize_t i;

	for (i = 0; i < n; i++, keys += elt_size) {
		hashset_remove(s, keys);
	}
}

const void *hashset_find(const struct hashset *s, const void *key,
			 struct hashset_pos *pos)
{
	assert(s);
	assert(key);
	assert(pos);
	assert(sparsetable_count(&s->table) < sparsetable_size(&s->table));

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

bool hashset_insert(struct hashset *s, struct hashset_pos *pos, const void *val)
{
	assert(!pos->has_existing);
	assert(hashset_hash(s, val) == pos->hash);
	assert(pos->has_insert);
	assert(hashset_size(s) < hashset_max_size(s));

	if (hashset_needs_grow_delta(s, 1)) {
		if (hashset_grow_delta(s, 1)) {
			hashset_find(s, val, pos);	// need to recompute pos
		} else {
			return false;
		}
	}

	return sparsetable_insert(&s->table, &pos->insert, val);
}

void hashset_replace(struct hashset *s, struct hashset_pos *pos,
		     const void *val)
{
	assert(pos->has_existing);
	assert(hashset_hash(s, val) == pos->hash);

	if (pos->has_insert && sparsetable_insert(&s->table, &pos->insert, val)) {
		sparsetable_erase(&s->table, &pos->existing);
	} else {
		// in case the insert fails, replace the existing value
		sparsetable_replace(&s->table, &pos->existing, val);
	}
}

void hashset_erase(struct hashset *s, struct hashset_pos *pos)
{
	assert(pos->has_existing);
	sparsetable_erase(&s->table, &pos->existing);
	if (hashset_needs_shrink(s)) {
		hashset_shrink(s);	// ok if this fails
	}
}

void hashset_iter_init(const struct hashset *s, struct hashset_iter *it)
{
	sparsetable_iter_init(&s->table, &it->table_it);
}

void hashset_iter_deinit(const struct hashset *s, struct hashset_iter *it)
{
	sparsetable_iter_deinit(&s->table, &it->table_it);
}

void hashset_iter_reset(const struct hashset *s, struct hashset_iter *it)
{
	sparsetable_iter_reset(&s->table, &it->table_it);
}

bool hashset_iter_advance(const struct hashset *s, struct hashset_iter *it)
{
	return sparsetable_iter_skip(&s->table, &it->table_it);
}

const void *hashset_iter_current(const struct hashset *s,
				 const struct hashset_iter *it)
{
	return sparsetable_iter_current(&s->table, &it->table_it);
}
