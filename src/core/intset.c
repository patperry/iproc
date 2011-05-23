#include "port.h"
#include <assert.h>
#include "intset.h"

DEFINE_COMPARE_FN(intptr_compare, intptr_t)

static intptr_t intset_at(const struct intset *s, ssize_t index);
static ssize_t intset_index(const struct intset *s, intptr_t val);

bool intset_init(struct intset *s)
{
	assert(s);

	if (array_init(&s->values, sizeof(intptr_t))) {
		return s;
	}
	return NULL;
}

bool intset_init_copy(struct intset *s, const struct intset *src)
{
	assert(s);
	assert(src);

	if (intset_init(s)) {
		if (intset_assign_copy(s, src)) {
			return s;
		}
		intset_deinit(s);
	}
	return NULL;
}

void intset_deinit(struct intset *s)
{
	assert(s);

	array_deinit(&s->values);
}

bool intset_assign_copy(struct intset *s, const struct intset *src)
{
	assert(s);
	assert(src);

	if (array_assign_copy(&s->values, &src->values)) {
		return s;
	}
	return NULL;
}

intptr_t *intset_copy_to(const struct intset *s, intptr_t *dst)
{
	assert(s);
	assert(dst || intset_empty(s));

	struct intset_iter it;
	intset_iter_init(s, &it);
	while (intset_iter_advance(s, &it)) {
		*dst++ = intset_iter_current(s, &it);
	}
	intset_iter_deinit(s, &it);
	return dst;
}

void intset_clear(struct intset *s)
{
	assert(s);
	array_clear(&s->values);
}

bool intset_empty(const struct intset *s)
{
	assert(s);
	return !intset_size(s);
}

ssize_t intset_size(const struct intset *s)
{
	assert(s);
	return array_count(&s->values);
}

intptr_t intset_min(const struct intset *s)
{
	assert(s);
	assert(!intset_empty(s));
	return *(intptr_t *)array_item(&s->values, 0);
}

intptr_t intset_max(const struct intset *s)
{
	assert(s);
	assert(!intset_empty(s));
	return *(intptr_t *)array_item(&s->values, array_count(&s->values) - 1);
}

bool intset_contains(const struct intset *s, intptr_t val)
{
	assert(s);

	struct intset_pos pos;
	return intset_find(s, val, &pos);
}

bool intset_add(struct intset *s, intptr_t val)
{
	assert(s);
	struct intset_pos pos;

	if (!intset_find(s, val, &pos)) {
		return intset_insert(s, &pos);
	}
	return true;		// already in set
}

bool intset_add_all(struct intset *s, const intptr_t *vals, ssize_t n)
{
	assert(s);
	assert(vals || n == 0);
	assert(n >= 0);

	ssize_t i;
	for (i = 0; i < n; i++) {
		if (!intset_add(s, vals[i]))
			goto rollback;
	}
	return true;
rollback:
	for (; i > 0; i--) {
		intset_remove(s, vals[i - 1]);
	}
	return false;
}

void intset_remove(struct intset *s, intptr_t val)
{
	assert(s);
	struct intset_pos pos;

	if (intset_find(s, val, &pos)) {
		intset_erase(s, &pos);
	}
}

void intset_remove_all(struct intset *s, const intptr_t *vals, ssize_t n)
{
	assert(s);
	assert(vals || n == 0);
	assert(n >= 0);

	ssize_t i;
	for (i = 0; i < n; i++) {
		intset_remove(s, vals[i]);
	}
}

bool intset_find(const struct intset *s, intptr_t val, struct intset_pos *pos)
{
	assert(s);
	assert(pos);

	ssize_t index = array_binary_search(&s->values, &val, intptr_compare);
	pos->value = val;

	if (index >= 0) {
		pos->index = index;
		return true;
	} else {
		pos->index = ~index;
		return false;
	}
}

bool intset_insert(struct intset *s, const struct intset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(pos->index == 0 || intset_at(s, pos->index - 1) < pos->value);
	assert(pos->index == intset_size(s)
	       || intset_at(s, pos->index) > pos->value);

	return array_insert(&s->values, pos->index, &pos->value);
}

void intset_erase(struct intset *s, const struct intset_pos *pos)
{
	assert(s);
	assert(pos);
	assert(0 <= pos->index && pos->index < intset_size(s));
	assert(intset_at(s, pos->index) == pos->value);

	array_remove_at(&s->values, pos->index);
}

void intset_iter_init(const struct intset *s, struct intset_iter *it)
{
	assert(s);
	assert(it);
	intset_iter_reset(s, it);
}

void intset_iter_deinit(const struct intset *s, struct intset_iter *it)
{
	assert(s);
	assert(it);
}

void intset_iter_reset(const struct intset *s, struct intset_iter *it)
{
	assert(s);
	assert(it);
	it->index = -1;
}

bool intset_iter_advance(const struct intset *s, struct intset_iter *it)
{
	assert(s);
	assert(it);
	it->index++;
	return it->index < intset_size(s);
}

intptr_t intset_iter_current(const struct intset *s,
			     const struct intset_iter *it)
{
	assert(s);
	assert(it);
	return intset_at(s, it->index);
}

intptr_t intset_at(const struct intset *s, ssize_t index)
{
	assert(s);
	assert(index >= 0);
	assert(index < intset_size(s));

	return *(intptr_t *)array_item(&s->values, index);
}

ssize_t intset_index(const struct intset *s, intptr_t val)
{
	assert(s);

	struct intset_pos pos;
	if (intset_find(s, val, &pos)) {
		return pos.index;
	}

	return -1;
}
