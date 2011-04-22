#include "port.h"
#include <assert.h>
#include "intmap.h"

static int intmap_item_compare(const void *x1, const void *x2);
static struct intmap_item *intmap_at(const struct intmap *m, ssize_t i);
static ssize_t intmap_index(const struct intmap *m, intptr_t key);


bool intmap_init(struct intmap *m)
{
	assert(m);
	return darray_init(&m->items, sizeof(struct intmap_item));
}

bool intmap_init_copy(struct intmap *m, const struct intmap * src)
{
	assert(m);
	assert(src);

	if (intmap_init(m)) {
		if (intmap_assign_copy(m, src)) {
			return true;
		}
		intmap_deinit(m);
	}
	return false;
}

void intmap_deinit(struct intmap *m)
{
	assert(m);

	darray_deinit(&m->items);
}

bool intmap_assign_copy(struct intmap *m, const struct intmap *src)
{
	assert(m);
	assert(src);
	return darray_assign_copy(&m->items, &src->items);
}

intptr_t *intmap_copy_vals_to(const struct intmap *m, intptr_t *dst)
{
	assert(m);
	assert(dst || intmap_empty(m));

	struct intmap_iter it;
	intmap_iter_init(m, &it);
	while (intmap_iter_advance(m, &it)) {
		*dst++ = intmap_iter_current(m, &it).val;
	}
	intmap_iter_deinit(m, &it);
	return dst;
}

intptr_t *intmap_copy_keys_to(const struct intmap * m, intptr_t * dst)
{
	assert(m);
	assert(dst || intmap_empty(m));
	
	struct intmap_iter it;
	intmap_iter_init(m, &it);
	while (intmap_iter_advance(m, &it)) {
		*dst++ = intmap_iter_current(m, &it).key;
	}
	intmap_iter_deinit(m, &it);
	return dst;
}

void intmap_clear(struct intmap *m)
{
	assert(m);
	darray_clear(&m->items);
}

bool intmap_empty(const struct intmap *m)
{
	assert(m);
	return darray_empty(&m->items);
}


ssize_t intmap_size(const struct intmap *m)
{
	assert(m);
	return darray_size(&m->items);
}

ssize_t intmap_max_size(const struct intmap *m)
{
	assert(m);
	return darray_max_size(&m->items);
}

bool intmap_contains(const struct intmap * m, intptr_t key)
{
	assert(m);
	struct intmap_pos pos;
	return intmap_find(m, key, &pos);
}

intptr_t intmap_lookup(const struct intmap *m, intptr_t key)
{
	assert(m);
	return intmap_lookup_with(m, key, (intptr_t)NULL);
}

intptr_t intmap_lookup_with(const struct intmap *m, intptr_t key,
			     intptr_t val0)
{
	assert(m);
	struct intmap_pos pos;
	const intptr_t *valp;
	
	if ((valp = intmap_find(m, key, &pos))) {
		return *valp;
	}

	return val0;
}

bool intmap_add(struct intmap *m, intptr_t key, intptr_t val)
{
	assert(m);
	assert(val);
	
	return intmap_add_all(m, &key, &val, 1);
}

bool intmap_add_all(struct intmap *m, const intptr_t *keys,
		    const intptr_t *vals, ssize_t n)
{
	assert(m);
	assert(keys || n == 0);
	assert(vals || n == 0);
	assert(n >= 0);

	struct intmap_pos pos;
	intptr_t oldvals[n];
	const intptr_t *valp;
	bool finds[n];
	ssize_t i;

	for (i = 0; i < n; i++) {
		valp = intmap_find(m, keys[i], &pos);
		finds[i] = valp;
		if ((finds[i])) {
			oldvals[i] = *valp;
			intmap_replace(m, &pos, vals[i]);
		} else {
			if (!intmap_insert(m, &pos, vals[i]))
				goto rollback;
		}
	}
	return true;
rollback:
	for (; i > 0; i--) {
		if (!finds[i-1]) {
			intmap_remove(m, keys[i-1]);
		} else {
			intmap_add(m, keys[i-1], oldvals[i-1]);
		}
	}
	return false;
}

void intmap_remove(struct intmap *m, intptr_t key)
{
	assert(m);
	intmap_remove_all(m, &key, 1);
}

void intmap_remove_all(struct intmap *m, const intptr_t * keys, ssize_t n)
{
	assert(m);
	assert(keys || n == 0);
	assert(n >= 0);

	struct intmap_pos pos;
	ssize_t i;

	for (i = 0; i < n; i++) {
		if (intmap_find(m, keys[i], &pos)) {
			darray_erase(&m->items, pos.index);
		}
	}
}

const intptr_t *intmap_find(const struct intmap *m, intptr_t key,
			    struct intmap_pos *pos)
{
	assert(m);
	assert(pos);

	ssize_t index = intmap_index(m, key);
	pos->key = key;
	
	if (index < 0) {
		pos->index = ~index;
		return NULL;
	} else {
		pos->index = index;
		return &intmap_at(m, pos->index)->val;
	}
}

bool intmap_insert(struct intmap *m, struct intmap_pos *pos,
		   intptr_t val)
{
	assert(m);
	assert(pos);
	assert(val);

	struct intmap_item item = { .key = pos->key, .val = val };

	return darray_insert(&m->items, pos->index, &item);
}

void intmap_replace(struct intmap *m, struct intmap_pos *pos, intptr_t val)
{
	assert(m);
	assert(pos);

	((struct intmap_item *)darray_at(&m->items, pos->index))->val = val;
}

void intmap_erase(struct intmap *m, struct intmap_pos *pos)
{
	assert(m);
	assert(pos);
	assert(((struct intmap_item *)darray_at(&m->items, pos->index))->key == pos->key);

	darray_erase(&m->items, pos->index);
}


void intmap_iter_init(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
	intmap_iter_reset(m, it);
}

void intmap_iter_deinit(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
}

void intmap_iter_reset(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
	it->index = -1;
}

bool intmap_iter_advance(const struct intmap *m, struct intmap_iter *it)
{
	assert(m);
	assert(it);
	
	it->index++;
	return (it->index < intmap_size(m));
}

struct intmap_item intmap_iter_current(const struct intmap *m,
				       const struct intmap_iter *it)
{
	assert(m);
	assert(it);
	
	return *intmap_at(m, it->index);
}

static int intmap_item_compare(const void *x1, const void *x2)
{
	assert(x1);
	assert(x2);
	
	const struct intmap_item *kv1 = x1;
	const struct intmap_item *kv2 = x2;

	return intptr_compare(&kv1->key, &kv2->key);
}

static struct intmap_item *intmap_at(const const struct intmap *m, ssize_t i)
{
	assert(m);
	assert(0 <= i && i < intmap_size(m));
	return darray_at(&m->items, i);
}

ssize_t intmap_index(const struct intmap *m, intptr_t key)
{
	assert(m);
	struct intmap_item kvkey = { .key = key };
	return darray_binary_search(&m->items, &kvkey, intmap_item_compare);
}
