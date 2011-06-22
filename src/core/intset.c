#include "port.h"
#include <assert.h>
#include "intset.h"

DEFINE_COMPARE_FN(intptr_compare, intptr_t)

static intptr_t intset_at(const struct intset *s, ssize_t index);

void intset_init(struct intset *s)
{
	assert(s);

	array_init(&s->keys, sizeof(intptr_t));
}

void intset_init_copy(struct intset *s, const struct intset *src)
{
	assert(s);
	assert(src);

	intset_init(s);
	intset_assign_copy(s, src);
}

void intset_assign_copy(struct intset *s, const struct intset *src)
{
	assert(s);
	assert(src);

	array_assign_copy(&s->keys, &src->keys);
}

void intset_deinit(struct intset *s)
{
	assert(s);

	array_deinit(&s->keys);
}

bool intset_add(struct intset *s, intptr_t key)
{
	assert(s);
	struct intset_pos pos;

	if (!intset_find(s, key, &pos)) {
		intset_insert(s, &pos);
		return true;
	}
	return false;	// already in set
}

void intset_clear(struct intset *s)
{
	assert(s);
	array_clear(&s->keys);
}

bool intset_contains(const struct intset *s, intptr_t key)
{
	assert(s);

	struct intset_pos pos;
	return intset_find(s, key, &pos);
}

void intset_copy_to(const struct intset *s, intptr_t *dst)
{
	assert(s);
	assert(dst || !intset_count(s));

	struct intset_iter it;
	INTSET_FOREACH(it, s) {
		*dst++ = INTSET_KEY(it);
	}
}

bool intset_remove(struct intset *s, intptr_t key)
{
	assert(s);
	bool found;
	struct intset_pos pos;

	if ((found = intset_find(s, key, &pos))) {
		intset_remove_at(s, &pos);
	}
	return found;
}

bool intset_find(const struct intset *s, intptr_t key, struct intset_pos *pos)
{
	assert(s);
	assert(pos);

	ssize_t index = array_binary_search(&s->keys, &key, intptr_compare);
	pos->key = key;

	if (index >= 0) {
		pos->index = index;
		return true;
	} else {
		pos->index = ~index;
		return false;
	}
}

void intset_insert(struct intset *s, const struct intset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(pos->index == 0 || intset_at(s, pos->index - 1) < pos->key);
	assert(pos->index == intset_count(s)
	       || intset_at(s, pos->index) > pos->key);

	array_insert(&s->keys, pos->index, &pos->key);
}

void intset_remove_at(struct intset *s, const struct intset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(0 <= pos->index && pos->index < intset_count(s));
	assert(intset_at(s, pos->index) == pos->key);

	array_remove_at(&s->keys, pos->index);
}

struct intset_iter intset_iter_make(const struct intset *s)
{
	assert(s);

	struct intset_iter it;
	it.set = s;
	intset_iter_reset(&it);
	return it;
}

void intset_iter_reset(struct intset_iter *it)
{
	assert(it);
	it->index = -1;
}

bool intset_iter_advance(struct intset_iter *it)
{
	assert(it);
	it->index++;
	return it->index < intset_count(it->set);
}

intptr_t intset_at(const struct intset *s, ssize_t index)
{
	assert(s);
	assert(index >= 0);
	assert(index < intset_count(s));

	return *(intptr_t *)array_item(&s->keys, index);
}
