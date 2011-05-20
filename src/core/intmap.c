#include "port.h"
#include <assert.h>
#include "intmap.h"

DEFINE_HASH_FN(intptr_hash, intptr_t)DEFINE_EQUALS_FN(intptr_equals, intptr_t)

bool intmap_init(struct intmap *m, size_t elt_size, size_t elt_align)
{
	assert(m);
	assert(elt_size >= 0);
	assert(elt_align > 0);

	size_t key_size = sizeof(intptr_t);
	size_t key_align = alignof(intptr_t);
	size_t val_offset =
	    ((char *)ALIGN_PTR((char *)NULL + key_size, elt_align) -
	     (char *)NULL);
	size_t pair_size =
	    ((char *)ALIGN_PTR((char *)NULL + val_offset + elt_size, key_align)
	     - (char *)NULL);

	m->elt_size = elt_size;
	m->elt_align = elt_align;
	m->val_offset = val_offset;

	return hashset_init(&m->pairs, intptr_hash, NULL, intptr_equals, NULL,
			    pair_size);
}

bool intmap_init_copy(struct intmap *m, const struct intmap *src)
{
	assert(m);
	assert(src);

	m->elt_size = src->elt_size;
	m->elt_align = src->elt_align;
	m->val_offset = src->val_offset;
	return hashset_init_copy(&m->pairs, &src->pairs);
}

void intmap_deinit(struct intmap *m)
{
	assert(m);
	hashset_deinit(&m->pairs);
}

bool intmap_assign_copy(struct intmap *m, const struct intmap *src)
{
	assert(m);
	assert(src);

	if (hashset_assign_copy(&m->pairs, &src->pairs)) {
		m->elt_size = src->elt_size;
		m->elt_align = src->elt_align;
		m->val_offset = src->val_offset;
		return true;

	}
	return false;
}

void *intmap_copy_vals_to(const struct intmap *m, void *dst)
{
	assert(m);
	assert(dst || intmap_empty(m));

	size_t elt_size = intmap_elt_size(m);
	struct intmap_iter it;
	const void *val;

	intmap_iter_init(m, &it);
	while (intmap_iter_advance(m, &it)) {
		val = intmap_iter_current(m, &it);
		memcpy(dst, val, elt_size);
		dst = (char *)dst + elt_size;
	}
	intmap_iter_deinit(m, &it);
	return dst;
}

intptr_t *intmap_copy_keys_to(const struct intmap *m, intptr_t *dst)
{
	assert(m);
	assert(dst || intmap_empty(m));

	struct intmap_iter it;
	intmap_iter_init(m, &it);
	while (intmap_iter_advance(m, &it)) {
		*dst++ = intmap_iter_current_key(m, &it);
	}
	intmap_iter_deinit(m, &it);
	return dst;
}

void intmap_clear(struct intmap *m)
{
	assert(m);
	hashset_clear(&m->pairs);
}

bool intmap_empty(const struct intmap *m)
{
	assert(m);
	return hashset_empty(&m->pairs);
}

ssize_t intmap_size(const struct intmap *m)
{
	assert(m);
	return hashset_size(&m->pairs);
}

size_t intmap_elt_size(const struct intmap *m)
{
	assert(m);
	return m->elt_size;
}

size_t intmap_elt_align(const struct intmap *m)
{
	assert(m);
	return m->elt_align;
}

bool intmap_contains(const struct intmap *m, intptr_t key)
{
	assert(m);
	struct intmap_pos pos;
	return intmap_find(m, key, &pos);
}

void *intmap_lookup(const struct intmap *m, intptr_t key)
{
	assert(m);
	return intmap_lookup_with(m, key, NULL);
}

void *intmap_lookup_with(const struct intmap *m, intptr_t key, const void *val0)
{
	assert(m);
	struct intmap_pos pos;
	void *val = intmap_find(m, key, &pos);
	return val ? val : (void *)val0;
}

bool intmap_add(struct intmap *m, intptr_t key, const void *val)
{
	assert(m);
	assert(val || intmap_elt_size(m) == 0);

	struct intmap_pos pos;
	if (intmap_find(m, key, &pos)) {
		intmap_replace(m, &pos, val);
		return true;
	} else {
		return intmap_insert(m, &pos, val);
	}
}

ssize_t intmap_add_all(struct intmap *m, const intptr_t *keys,
		       const void *vals, ssize_t n)
{
	assert(m);
	assert(keys || n == 0);
	assert(vals || n == 0 || intmap_elt_size(m) == 0);
	assert(n >= 0);

	size_t elt_size = intmap_elt_size(m);
	ssize_t i;

	for (i = 0; i < n; i++) {
		if (!intmap_add(m, *keys, vals))
			break;
		keys++;
		vals = (char *)vals + elt_size;
	}
	return i;
}

void intmap_remove(struct intmap *m, intptr_t key)
{
	assert(m);

	struct intmap_pos pos;
	if (intmap_find(m, key, &pos)) {
		intmap_erase(m, &pos);
	}
}

void intmap_remove_all(struct intmap *m, const intptr_t *keys, ssize_t n)
{
	assert(m);
	assert(keys || n == 0);
	assert(n >= 0);

	ssize_t i;

	for (i = 0; i < n; i++) {
		intmap_remove(m, keys[i]);
	}
}

void *intmap_find(const struct intmap *m, intptr_t key, struct intmap_pos *pos)
{
	assert(m);
	assert(pos);

	void *pair;

	pos->key = key;
	if ((pair = hashset_find(&m->pairs, &key, &pos->pairs_pos))) {
		return (char *)pair + m->val_offset;
	}
	return NULL;
}

void *intmap_insert(struct intmap *m, struct intmap_pos *pos, const void *val)
{
	assert(m);
	assert(pos);

	void *pair = hashset_insert(&m->pairs, &pos->pairs_pos, &pos->key);
	void *res = NULL;

	if (pair) {
		res = (char *)pair + m->val_offset;
		if (val)
			memcpy(res, val, m->elt_size);
	}
	return res;
}

void *intmap_replace(struct intmap *m, struct intmap_pos *pos, const void *val)
{
	assert(m);
	assert(pos);

	void *pair = hashset_replace(&m->pairs, &pos->pairs_pos, &pos->key);
	void *res = NULL;

	if (pair) {
		res = (char *)pair + m->val_offset;
		if (val)
			memcpy(res, val, m->elt_size);
	}
	return res;
}

void intmap_erase(struct intmap *m, struct intmap_pos *pos)
{
	assert(m);
	assert(pos);

	hashset_erase(&m->pairs, &pos->pairs_pos);
}

void intmap_iter_init(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
	hashset_iter_init(&m->pairs, &it->pairs_it);
}

void intmap_iter_deinit(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
	hashset_iter_deinit(&m->pairs, &it->pairs_it);
}

void intmap_iter_reset(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
	hashset_iter_reset(&m->pairs, &it->pairs_it);
}

bool intmap_iter_advance(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);

	return hashset_iter_advance(&m->pairs, &it->pairs_it);
}

void *intmap_iter_current(const struct intmap *m, const struct intmap_iter *it)
{
	assert(m);
	assert(it);

	void *pair = hashset_iter_current(&m->pairs, &it->pairs_it);
	return (char *)pair + m->val_offset;
}

intptr_t intmap_iter_current_key(const struct intmap *m,
				 const struct intmap_iter *it)
{
	assert(m);
	assert(it);

	void *pair = hashset_iter_current(&m->pairs, &it->pairs_it);
	return *(intptr_t *)pair;
}
