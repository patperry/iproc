#ifndef HASHSET_H
#define HASHSET_H

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

#include <stddef.h>

struct hashset {
	size_t width;
	size_t (*hash) (const void *);
	int (*compar) (const void *, const void *);

	size_t nbucket;
	void *buckets;
	unsigned char *status;

	size_t count;
	size_t count_max;
};

struct hashset_pos {
	size_t insert;
	size_t existing;
};

struct hashset_iter {
	const struct hashset *s;
	size_t i;
	void *val;
};

#define HASHSET_VAL(it) ((it).val)
#define HASHSET_FOREACH(it, set) \
	for ((it) = hashset_iter_make(set); hashset_iter_advance(&(it));)

/* create, destroy */
void hashset_init(struct hashset *s, size_t width,
		  size_t (*hash) (const void *),
		  int (*compar) (const void *, const void *));
void hashset_init_copy(struct hashset *s, const struct hashset *src);
void hashset_assign_copy(struct hashset *s, const struct hashset *src);
void hashset_deinit(struct hashset *s);

/* properties */
static inline size_t hashset_count(const struct hashset *s);
static inline size_t hashset_width(const struct hashset *s);
static inline int hashset_compare(const struct hashset *s, const void *key1,
				  const void *key2);
static inline size_t hashset_hash(const struct hashset *s, const void *key);

static inline size_t hashset_capacity(const struct hashset *s);
void hashset_ensure_capacity(struct hashset *s, size_t n);

void *hashset_item(const struct hashset *s, const void *key);
void *hashset_set_item(struct hashset *s, const void *val);

/* methods */
void *hashset_add(struct hashset *s, const void *val);
void hashset_clear(struct hashset *s);
int hashset_contains(const struct hashset *s, const void *key);
int hashset_remove(struct hashset *s, const void *key);
void hashset_trim_excess(struct hashset *s);

/* position-based operations */
void *hashset_find(const struct hashset *s, const void *key,
		   struct hashset_pos *pos);
void *hashset_insert(struct hashset *s, struct hashset_pos *pos,
		     const void *val);
void hashset_remove_at(struct hashset *s, struct hashset_pos *pos);

/* iteration */
struct hashset_iter hashset_iter_make(const struct hashset *s);
void hashset_iter_reset(struct hashset_iter *it);
void *hashset_iter_advance(struct hashset_iter *it);

/* static method definitions */
size_t hashset_count(const struct hashset *s)
{
	return s->count;
}

static inline size_t hashset_capacity(const struct hashset *s)
{
	return s->count_max;
}

size_t hashset_width(const struct hashset *s)
{
	return s->width;
}

static inline int hashset_compare(const struct hashset *s, const void *key1,
				  const void *key2)
{
	return s->compar(key1, key2);
}

static inline size_t hashset_hash(const struct hashset *s, const void *key)
{
	return s->hash(key);
}

#endif /* HASHSET_H */
