#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "util.h"
#include "shashset.h"

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
static void shashset_reset_thresholds(struct shashset *s, ssize_t num_buckets)
{
	s->enlarge_threshold = num_buckets * HT_ENLARGE_FACTOR;
}

static ssize_t shashset_bucket_count(const struct shashset *s)
{
	return sparsetable_size(&s->table);
}

static void shashset_init_sized(struct shashset *s, hash_fn hash,
			       equals_fn equals, ssize_t num_buckets,
			       size_t elt_size)
{
	assert(s);
	assert(hash);
	assert(equals);
	assert(num_buckets >= HT_MIN_BUCKETS);

	sparsetable_init(&s->table, num_buckets, elt_size);
	s->hash = hash;
	s->equals = equals;
	shashset_reset_thresholds(s, num_buckets);
}

static void shashset_init_copy_sized(struct shashset *s,
				    const struct shashset *src,
				    ssize_t num_buckets)
{
	assert(s);
	assert(src);
	assert(s != src);
	assert(num_buckets >= HT_MIN_BUCKETS);

	struct shashset_iter it;
	const void *key;

	shashset_init_sized(s, src->hash, src->equals, num_buckets,
			   shashset_elt_size(src));

	SHASHSET_FOREACH(it, src) {
		key = SHASHSET_KEY(it);
		shashset_add(s, key);
	}
}

static bool shashset_needs_grow_delta(const struct shashset *s, ssize_t delta)
{
	assert(delta >= 0);
	assert(sparsetable_count(&s->table) <= HT_MAX_SIZE - delta);

	if (shashset_bucket_count(s) >= HT_MIN_BUCKETS
	    && sparsetable_count(&s->table) + delta <= s->enlarge_threshold) {
		return false;
	} else {
		return true;
	}
}

/* the implementation is simpler than in sparsehashtable.h since we don't
 * store deleted keys */
static void shashset_grow_delta(struct shashset *s, ssize_t delta)
{
	ssize_t num_nonempty = sparsetable_count(&s->table);
	ssize_t bucket_count = shashset_bucket_count(s);
	ssize_t resize_to = min_buckets(num_nonempty + delta, bucket_count);

	struct shashset copy;
	shashset_init_copy_sized(&copy, s, resize_to);
	shashset_deinit(s);
	*s = copy;
}

void shashset_init(struct shashset *s, hash_fn hash, equals_fn equals,
		  size_t elt_size)
{
	assert(s);
	assert(hash);
	assert(equals);

	ssize_t num_buckets = HT_DEFAULT_STARTING_BUCKETS;
	shashset_init_sized(s, hash, equals, num_buckets, elt_size);
}

void shashset_init_copy(struct shashset *s, const struct shashset *src)
{
	assert(s);
	assert(src);
	assert(s != src);

	ssize_t num_buckets = shashset_bucket_count(src);
	shashset_init_copy_sized(s, src, num_buckets);
}

void shashset_assign_copy(struct shashset *s, const struct shashset *src)
{
	assert(s);
	assert(src);

	shashset_deinit(s);
	shashset_init_copy(s, src);
}

void shashset_deinit(struct shashset *s)
{
	assert(s);

	sparsetable_deinit(&s->table);
}

void *shashset_item(const struct shashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct shashset_pos pos;
	return shashset_find(s, key, &pos);
}

void *shashset_set_item(struct shashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct shashset_pos pos;
	void *dst;

	if ((dst = shashset_find(s, key, &pos))) {
		memcpy(dst, key, shashset_elt_size(s));
		return dst;
	} else {
		return shashset_insert(s, &pos, key);
	}
}

void *shashset_add(struct shashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct shashset_pos pos;
	if (shashset_find(s, key, &pos)) {
		return NULL;
	} else {
		return shashset_insert(s, &pos, key);
	}
}

void shashset_clear(struct shashset *s)
{
	assert(s);

	ssize_t num_buckets = HT_DEFAULT_STARTING_BUCKETS;

	sparsetable_clear(&s->table);
	sparsetable_set_size(&s->table, num_buckets);
	shashset_reset_thresholds(s, num_buckets);
}

bool shashset_contains(const struct shashset *s, const void *key)
{
	assert(s);
	assert(key);

	return shashset_item(s, key);
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

bool shashset_remove(struct shashset *s, const void *key)
{
	assert(s);
	assert(key);

	bool found;

	struct shashset_pos pos;
	if ((found = shashset_find(s, key, &pos))) {
		shashset_remove_at(s, &pos);
	}
	assert(!shashset_contains(s, key));
	return found;
}

/* MISSING remove_where */
/* MISSING set_equals */
/* MISSING symmetric_except_with */
/* MISSING to_string */

void shashset_trim_excess(struct shashset *s)
{
	assert(s);

	ssize_t count = sparsetable_count(&s->table);
	ssize_t resize_to = min_buckets(count, 0);

	struct shashset copy;
	shashset_init_copy_sized(&copy, s, resize_to);
	shashset_deinit(s);
	*s = copy;
}

/* MISSING union_with */

void *shashset_find(const struct shashset *s, const void *key,
		   struct shashset_pos *pos)
{
	assert(s);
	assert(key);
	assert(pos);

	const struct sparsetable *table = &s->table;
	const ssize_t bucket_count = sparsetable_size(table);
	ssize_t num_probes = 0;	// how many times we've probed
	const ssize_t bucket_count_minus_one = bucket_count - 1;
	uint32_t hash = shashset_hash(s, key);
	ssize_t bucknum = hash & bucket_count_minus_one;
	struct sparsetable_pos table_pos;
	void *ptr;
	bool deleted;

	pos->hash = hash;
	pos->has_insert = false;
	pos->has_existing = false;

	for (num_probes = 0; num_probes < bucket_count; num_probes++) {
		ptr = sparsetable_find(table, bucknum, &table_pos);
		deleted = sparsetable_deleted(table, &table_pos);
		if (!ptr && !deleted) {	// bucket is empty
			if (!pos->has_insert) {	// found no prior place to insert
				pos->insert = table_pos;
				pos->has_insert = true;
			}
			pos->has_existing = false;
			return NULL;
		} else if (!ptr && deleted) {	// keep searching, but mark to insert
			if (!pos->has_insert) {
				pos->insert = table_pos;
				pos->has_insert = true;
			}
		} else if (shashset_equals(s, key, ptr)) {
			pos->existing = table_pos;
			pos->has_existing = true;
			return ptr;
		}
		bucknum =
		    (bucknum +
		     JUMP_(key, num_probes + 1)) & bucket_count_minus_one;
	}

	return NULL;		// table is full and key is not present
}

void *shashset_insert(struct shashset *s, struct shashset_pos *pos,
		     const void *key)
{
	assert(s);
	assert(pos);
	assert(!pos->has_existing);
	assert(key);
	assert(shashset_hash(s, key) == pos->hash);

	if (shashset_needs_grow_delta(s, 1)) {
		shashset_grow_delta(s, 1);
		shashset_find(s, key, pos);	// need to recompute pos
	}

	assert(pos->has_insert);
	return sparsetable_insert(&s->table, &pos->insert, key);
}

void shashset_remove_at(struct shashset *s, struct shashset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(pos->has_existing);

	sparsetable_remove_at(&s->table, &pos->existing);
}

struct shashset_iter shashset_iter_make(const struct shashset *s)
{
	assert(s);

	struct shashset_iter it;
	it.table_it = sparsetable_iter_make(&s->table);
	return it;
}

void shashset_iter_reset(struct shashset_iter *it)
{
	assert(it);

	sparsetable_iter_reset(&it->table_it);
}

bool shashset_iter_advance(struct shashset_iter *it)
{
	assert(it);

	return sparsetable_iter_advance(&it->table_it);
}
