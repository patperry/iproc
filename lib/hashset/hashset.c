/* Copyright (c) 2011, Patrick O. Perry
 * Copyright (c) 2005, Google Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following disclaimer
 *       in the documentation and/or other materials provided with the
 *       distribution.
 *     * Neither the names of Patrick O. Perry, Google Inc,. nor the names of
 *       their contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>		// assert
#include <limits.h>		// CHAR_BIT
#include <stddef.h>		// size_t, NULL
#include <stdlib.h>		// free
#include <string.h>		// memset, memcpy
#include "xalloc.h"		// xrealloc

#include "hashset.h"

#define  HT_BUCKET_FULL    1
#define  HT_BUCKET_DELETED 2

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

/* The number of buckets must be a power of 2.  This is the largest 
 * power of 2 that a size_t can hold.
 */
#define HT_MAX_BUCKETS	((size_t)1 << (CHAR_BIT * sizeof(size_t) - 1))

/* The following is more accurate than ((pct)/100.0 * (x)) when x is
 * really big (> 2^52).
 */
#define PERCENT(pct,x) \
	((pct) * ((x) / 100) \
	 + (size_t)((pct) * (((x) % 100) / 100.0)))

#define HT_MAX_COUNT	PERCENT(HT_OCCUPANCY_PCT, HT_MAX_BUCKETS)

/* This is the smallest size a hashtable can be without being too crowded
 * If you like, you can give a min #buckets as well as a min #elts */
static size_t min_buckets(size_t count, size_t nbucket0)
{
	assert(count <= HT_MAX_COUNT);
	assert(nbucket0 <= HT_MAX_BUCKETS);

	size_t n = HT_MIN_BUCKETS;	// min buckets allowed

	while (n < nbucket0 || count > PERCENT(HT_OCCUPANCY_PCT, n)) {
		assert(2 * n > n);
		n *= 2;
	}

	assert(n >= nbucket0);
	assert(count <= PERCENT(HT_OCCUPANCY_PCT, n));

	return n;
}

/* Reset the enlarge threshold */
static void hashset_reset_thresholds(struct hashset *s, size_t nbucket)
{
	s->count_max = PERCENT(HT_OCCUPANCY_PCT, nbucket);
}

static size_t hashset_bucket_count(const struct hashset *s)
{
	return s->nbucket;
}

static void hashset_init_sized(struct hashset *s,
			       size_t width,
			       size_t (*hash) (const void *),
			       int (*compar) (const void *, const void *),
			       size_t nbucket)
{
	assert(s);
	assert(hash);
	assert(compar);
	assert(nbucket >= HT_MIN_BUCKETS);

	hashset_init(s, width, hash, compar);

	s->nbucket = nbucket;
	s->buckets = xcalloc(nbucket, width);
	s->status = xcalloc(nbucket, sizeof(s->status[0]));

	hashset_reset_thresholds(s, nbucket);
}

static void hashset_init_copy_sized(struct hashset *s,
				    const struct hashset *src, size_t nbucket)
{
	assert(s);
	assert(src);
	assert(s != src);
	assert(nbucket >= HT_MIN_BUCKETS);

	struct hashset_iter it;
	const void *val;

	hashset_init_sized(s, src->width, src->hash, src->compar, nbucket);

	HASHSET_FOREACH(it, src) {
		val = HASHSET_VAL(it);
		hashset_add(s, val);
	}
}

static int hashset_needs_grow_delta(const struct hashset *s, size_t delta)
{
	assert(delta <= HT_MAX_COUNT);
	assert(s->nbucket <= HT_MAX_COUNT - delta);
	assert(s->count_max >= s->count);

	if (s->nbucket >= HT_MIN_BUCKETS && delta <= s->count_max - s->count) {
		return 0;
	} else {
		return 1;
	}
}

static void hashset_grow_delta(struct hashset *s, size_t delta)
{
	assert(delta <= HT_MAX_COUNT - s->nbucket);

	size_t count0 = s->count;
	size_t nbucket0 = s->nbucket;
	size_t count = count0 + delta;
	size_t nbucket = min_buckets(count, nbucket0);

	if (nbucket > nbucket0) {
		struct hashset copy;
		hashset_init_copy_sized(&copy, s, nbucket);
		hashset_deinit(s);
		*s = copy;
	}
}

void hashset_ensure_capacity(struct hashset *s, size_t n)
{
	assert(s);
	assert(n >= hashset_count(s));
	assert(hashset_capacity(s) >= hashset_count(s));
	assert(n <= HT_MAX_BUCKETS);

	if (n > hashset_capacity(s)) {
		size_t delta = n - hashset_count(s);
		hashset_grow_delta(s, delta);
	}
}

void hashset_init(struct hashset *s,
		  size_t width,
		  size_t (*hash) (const void *),
		  int (*compar) (const void *, const void *))
{
	assert(s);
	assert(hash);
	assert(compar);

	s->buckets = NULL;
	s->nbucket = 0;
	s->width = width;
	s->status = NULL;
	s->count = 0;
	s->hash = hash;
	s->compar = compar;
	hashset_reset_thresholds(s, 0);
}

void hashset_init_copy(struct hashset *s, const struct hashset *src)
{
	assert(s);
	assert(src);
	assert(s != src);

	size_t nbucket = hashset_bucket_count(src);
	hashset_init_copy_sized(s, src, nbucket);
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

	free(s->status);
	free(s->buckets);
}

void *hashset_item(const struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	const void *buckets = s->buckets;
	const unsigned char *status = s->status;
	const size_t bucket_count = s->nbucket;
	const size_t width = s->width;
	size_t num_probes = 0;	// how many times we've probed
	const size_t bucket_count_minus_one = bucket_count - 1;
	size_t hash = hashset_hash(s, key);
	size_t bucknum = hash & bucket_count_minus_one;
	void *ptr;
	unsigned char stat;

	for (num_probes = 0; num_probes < bucket_count; num_probes++) {
		ptr = (char *)buckets + bucknum * width;
		stat = status[bucknum];
		if (!stat) {	// bucket is empty
			return NULL;
		} else if (stat == HT_BUCKET_DELETED) {	// bucket empty, deleted; keep searching

		} else if (!hashset_compare(s, key, ptr)) {
			return ptr;
		}
		bucknum =
		    (bucknum +
		     JUMP_(key, num_probes + 1)) & bucket_count_minus_one;
	}

	return NULL;
}

void *hashset_set_item(struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	struct hashset_pos pos;
	void *dst;

	if ((dst = hashset_find(s, key, &pos))) {
		return memcpy(dst, key, s->width);
	} else {
		return hashset_insert(s, &pos, key);
	}
}

void *hashset_add(struct hashset *s, const void *val)
{
	assert(s);
	assert(val);

	struct hashset_pos pos;
	if (hashset_find(s, val, &pos)) {
		return NULL;
	} else {
		return hashset_insert(s, &pos, val);
	}
}

void hashset_clear(struct hashset *s)
{
	assert(s);

	size_t n = hashset_bucket_count(s);

	memset(s->buckets, 0, n * s->width);
	memset(s->status, 0, n * sizeof(s->status[0]));
	s->count = 0;
}

int hashset_contains(const struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	if (hashset_item(s, key)) {
		return 1;
	} else {
		return 0;
	}
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

int hashset_remove(struct hashset *s, const void *key)
{
	assert(s);
	assert(key);

	const void *buckets = s->buckets;
	unsigned char *status = s->status;
	const size_t bucket_count = s->nbucket;
	const size_t width = s->width;
	size_t num_probes = 0;	// how many times we've probed
	const size_t bucket_count_minus_one = bucket_count - 1;
	size_t hash = hashset_hash(s, key);
	size_t bucknum = hash & bucket_count_minus_one;
	void *ptr;
	unsigned char stat;

	for (num_probes = 0; num_probes < bucket_count; num_probes++) {
		ptr = (char *)buckets + bucknum * width;
		stat = status[bucknum];
		if (!stat) {	// bucket is empty
			return 0;
		} else if (stat == HT_BUCKET_DELETED) {	// empty, deleted; keep searching
		} else if (!hashset_compare(s, key, ptr)) {
			status[bucknum] = HT_BUCKET_DELETED;
			s->count--;
			return 1;
		}
		bucknum =
		    (bucknum +
		     JUMP_(key, num_probes + 1)) & bucket_count_minus_one;
	}

	return 0;
}

/* MISSING remove_where */
/* MISSING set_equals */
/* MISSING symmetric_except_with */
/* MISSING to_string */

void hashset_trim_excess(struct hashset *s)
{
	assert(s);

	size_t count = hashset_count(s);
	size_t resize_to = min_buckets(count, 0);

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

	const void *buckets = s->buckets;
	const unsigned char *status = s->status;
	const size_t bucket_count = s->nbucket;
	const size_t width = s->width;
	size_t num_probes = 0;	// how many times we've probed
	const size_t bucket_count_minus_one = bucket_count - 1;
	size_t hash = hashset_hash(s, key);
	size_t bucknum = hash & bucket_count_minus_one;
	void *ptr;
	unsigned char stat;

	pos->insert = HT_MAX_BUCKETS;
	pos->existing = HT_MAX_BUCKETS;

	for (num_probes = 0; num_probes < bucket_count; num_probes++) {
		ptr = (char *)buckets + bucknum * width;
		stat = status[bucknum];
		if (!stat) {	// bucket is empty
			if (pos->insert == HT_MAX_BUCKETS) {	// found no prior place to insert
				pos->insert = bucknum;
			}
			return NULL;
		} else if (stat == HT_BUCKET_DELETED) {	// empty, deleted; keep searching, mark to insert
			if (pos->insert == HT_MAX_BUCKETS) {
				pos->insert = bucknum;
			}
		} else if (!hashset_compare(s, key, ptr)) {
			pos->existing = bucknum;
			return ptr;
		}
		bucknum =
		    (bucknum +
		     JUMP_(key, num_probes + 1)) & bucket_count_minus_one;
	}

	return NULL;		// table is full and key is not present
}

void *hashset_insert(struct hashset *s, struct hashset_pos *pos,
		     const void *val)
{
	assert(s);
	assert(pos);
	assert(pos->existing == HT_MAX_BUCKETS);
	assert(val);

	if (hashset_needs_grow_delta(s, 1)) {
		hashset_grow_delta(s, 1);
		hashset_find(s, val, pos);	// need to recompute pos
	}

	assert(!hashset_needs_grow_delta(s, 1));
	assert(s->count < s->count_max);
	assert(pos->insert != HT_MAX_BUCKETS);

	size_t ix = pos->insert;
	size_t width = s->width;

	pos->existing = ix;
	s->count++;
	s->status[ix] |= HT_BUCKET_FULL;

	void *ptr = s->buckets + ix * width;
	memcpy(ptr, val, width);

	return ptr;
}

void hashset_remove_at(struct hashset *s, struct hashset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(pos->existing != HT_MAX_BUCKETS);

	size_t ix = pos->existing;

	pos->insert = ix;
	pos->existing = HT_MAX_BUCKETS;

	s->count--;
	s->status[ix] = HT_BUCKET_DELETED;
}

struct hashset_iter hashset_iter_make(const struct hashset *s)
{
	assert(s);

	struct hashset_iter it;
	it.s = s;
	hashset_iter_reset(&it);
	return it;
}

void hashset_iter_reset(struct hashset_iter *it)
{
	assert(it);

	it->i = 0;
	it->val = NULL;
}

void *hashset_iter_advance(struct hashset_iter *it)
{
	assert(it);

	const struct hashset *s = it->s;
	const unsigned char *status = it->s->status;
	size_t i, n = hashset_bucket_count(s);

	for (i = it->i; i < n; i++) {
		if (status[i] & HT_BUCKET_FULL) {
			it->val = s->buckets + i * s->width;
			goto out;
		}
	}
	it->val = NULL;
out:
	it->i = i + 1;
	return it->val;
}
