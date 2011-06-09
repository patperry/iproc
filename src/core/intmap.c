#include "port.h"
#include <assert.h>
#include "intmap.h"

DEFINE_HASH_FN(intptr_hash, intptr_t)
DEFINE_EQUALS_FN(intptr_equals, intptr_t)


void intmap_init(struct intmap *m, size_t elt_size, size_t elt_align)
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

	hashset_init(&m->pairs, intptr_hash, intptr_equals, pair_size);
}

void intmap_init_copy(struct intmap *m, const struct intmap *src)
{
	assert(m);
	assert(src);

	m->elt_size = src->elt_size;
	m->elt_align = src->elt_align;
	m->val_offset = src->val_offset;
	hashset_init_copy(&m->pairs, &src->pairs);
}

void intmap_assign_copy(struct intmap *m, const struct intmap *src)
{
	assert(m);
	assert(src);

	hashset_assign_copy(&m->pairs, &src->pairs);
	m->elt_size = src->elt_size;
	m->elt_align = src->elt_align;
	m->val_offset = src->val_offset;
}

void intmap_deinit(struct intmap *m)
{
	assert(m);
	hashset_deinit(&m->pairs);
}

void *intmap_item(const struct intmap *m, intptr_t key)
{
	assert(m);
	struct intmap_pos pos;
	return intmap_find(m, key, &pos);
}

void *intmap_set_item(struct intmap *m, intptr_t key, const void *val)
{
	assert(m);

	struct intmap_pos pos;
	void *dst;

	if ((dst = intmap_find(m, key, &pos))) {
		assert(*(intptr_t *)((char *)dst - m->val_offset) == key);
		if (val) {
			memcpy(dst, val, intmap_elt_size(m));
		} else {
			memset(dst, 0, intmap_elt_size(m));
		}
		return dst;
	} else {
		return intmap_insert(m, &pos, val);
	}
}

void *intmap_add(struct intmap *m, intptr_t key, const void *val)
{
	assert(m);
	assert(val || intmap_elt_size(m) == 0);
	assert(!intmap_contains_key(m, key));

	struct intmap_pos pos;

	intmap_find(m, key, &pos);
	return intmap_insert(m, &pos, val);
}

void intmap_clear(struct intmap *m)
{
	assert(m);
	hashset_clear(&m->pairs);
}

bool intmap_contains_key(const struct intmap *m, intptr_t key)
{
	assert(m);
	struct intmap_pos pos;
	return intmap_find(m, key, &pos);
}

bool intmap_contains_val(const struct intmap *m, predicate_fn match_val,
			 void *udata)
{
	assert(m);
	assert(match_val);

	struct intmap_iter it;
	INTMAP_FOREACH(it, m) {
		if (match_val(INTMAP_VAL(it), udata))
			return true;
	}

	return false;
}

void intmap_copy_keys_to(const struct intmap *m, intptr_t *dst)
{
	assert(m);
	assert(dst || !intmap_count(m));

	struct intmap_iter it;

	INTMAP_FOREACH(it, m) {
		*dst++ = INTMAP_KEY(it);
	}
}

void intmap_copy_vals_to(const struct intmap *m, void *dst)
{
	assert(m);
	assert(dst || !intmap_count(m));

	size_t elt_size = intmap_elt_size(m);
	struct intmap_iter it;
	const void *val;

	INTMAP_FOREACH(it, m) {
		val = INTMAP_VAL(it);
		memcpy(dst, val, elt_size);
		dst = (char *)dst + elt_size;
	}
}

bool intmap_remove(struct intmap *m, intptr_t key)
{
	assert(m);

	struct intmap_pos pos;
	bool found;

	if ((found = intmap_find(m, key, &pos))) {
		intmap_remove_at(m, &pos);
	}

	return found;
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
		if (val) {
			memcpy(res, val, m->elt_size);
		} else {
			memset(res, 0, m->elt_size);
		}
	}
	return res;
}

void intmap_remove_at(struct intmap *m, struct intmap_pos *pos)
{
	assert(m);
	assert(pos);

	hashset_remove_at(&m->pairs, &pos->pairs_pos);
}

struct intmap_iter intmap_iter_make(const struct intmap *m)
{
	assert(m);
	struct intmap_iter it;
	it.map = m;
	it.pairs_it = hashset_iter_make(&m->pairs);
	return it;
}

void intmap_iter_reset(struct intmap_iter *it)
{
	assert(it);
	hashset_iter_reset(&it->pairs_it);
}

bool intmap_iter_advance(struct intmap_iter *it)
{
	assert(it);
	return hashset_iter_advance(&it->pairs_it);
}
