#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "hashset.h"

/* The probing method */
/* #define JUMP_(key, num_probes) (1)          // Linear probing */
#define JUMP_(key, num_probes)    (num_probes) // Quadratic probing

/* How full we let the table get before we resize, by default.
 * Knuth says .8 is good -- higher causes us to probe too much,
 * though it saves memory.
 */
#define HT_OCCUPANCY_PCT 80 // (out of 100);

/* How empty we let the table get before we resize lower, by default.
 * (0.0 means never resize lower.)
 * It should be less than OCCUPANCY_PCT / 2 or we thrash resizing
 */
#define HT_EMPTY_PCT	(0.4f * HT_OCCUPANCY_PCT)

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

#define HT_ENLARGE_FACTOR	(HT_OCCUPANCY_PCT / 100.0f)
#define HT_SHRINK_FACTOR	(HT_EMPTY_PCT / 100.0f)

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
	
	float enlarge = HT_ENLARGE_FACTOR;
	ssize_t sz = HT_MIN_BUCKETS;             // min buckets allowed

	while (sz < min_buckets_wanted
	       || num_elts > (ssize_t)(sz * enlarge)) {
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

static bool hashset_shrink(const struct hashset *s)
{
	assert(hashset_needs_shrink(s));

	ssize_t num_remain = sparsetable_count(&s->table);
	ssize_t bucket_count = hashset_bucket_count(s);
	float shrink_factor = HT_SHRINK_FACTOR;
	ssize_t sz = bucket_count / 2;    // find how much we should shrink

	while (sz > HT_DEFAULT_STARTING_BUCKETS &&
	       num_remain < (ssize_t)(sz * shrink_factor)) {
		sz /= 2;		// stay a power of 2
	}

		// TODO: now do the shrink
	return true;
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
static bool hashset_grow_delta(const struct hashset *s, ssize_t delta)
{
	ssize_t num_nonempty = sparsetable_count(&s->table);
	ssize_t bucket_count = hashset_bucket_count(s);
	ssize_t resize_to = min_buckets(num_nonempty + delta, bucket_count);

		// TODO: now do the grow
	return resize_to;
}

void *hashset_find(const struct hashset *s, const void *key,
		   struct hashset_pos *pos)
{
	assert(s);
	assert(key);
	assert(pos);
	assert(sparsetable_count(&s->table) < sparsetable_size(&s->table));

	const struct sparsetable *table = &s->table;
	const ssize_t bucket_count = sparsetable_size(table);
	ssize_t num_probes = 0;              // how many times we've probed
	const ssize_t bucket_count_minus_one = bucket_count - 1;
	uint32_t hash = hashset_hash(s, key);
	ssize_t bucknum = hash & bucket_count_minus_one;
	struct sparsetable_pos table_pos;
	void *val;
	bool deleted;

	pos->hash = hash;
	pos->has_insert = false;
	pos->has_existing = false;
	
	while (true) {                          // probe until something happens
		val = sparsetable_find(table, bucknum, &table_pos);
		deleted = sparsetable_deleted(table, &table_pos);
		if (!val && !deleted) { // bucket is empty
			if (!pos->has_insert) { // found no prior place to insert
				pos->insert = table_pos;
				pos->has_insert = true;
			}
			pos->has_existing = false;
			return NULL;
		} else if (!val && deleted) {// keep searching, but mark to insert
			if (!pos->has_insert) {
				pos->insert = table_pos;
				pos->has_insert = true;
			}
		} else if (hashset_equals(s, key, val)) {
			pos->existing = table_pos;
			pos->has_existing = true;
			return val;
		}
		num_probes++;                        // we're doing another probe
		bucknum = (bucknum + JUMP_(key, num_probes)) & bucket_count_minus_one;
		assert(num_probes < bucket_count
		       && "Hashtable is full: an error in hash or equals");
	}
}

bool hashset_insert(struct hashset *s, struct hashset_pos *pos,
		    const void *val)
{
	assert(!pos->has_existing);
	assert(hashset_hash(s, val) == pos->hash);
	assert(pos->has_insert);
	assert(hashset_size(s) < HT_MAX_SIZE);

	if (hashset_needs_grow_delta(s, 1)) {
		if (hashset_grow_delta(s, 1)) {
			hashset_find(s, val, pos); // need to recompute pos
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
	
	if (pos->has_insert
	    && sparsetable_insert(&s->table, &pos->insert, val)) {
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
		hashset_shrink(s); // ok if this fails
	}
}
