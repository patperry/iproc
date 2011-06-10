#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "util.h"
#include "cohort.h"

struct cohort *cohort_alloc(const struct vector *x)
{
	assert(x);

	struct cohort *c = xcalloc(1, sizeof(*c));
	cohort_init(c, x);
	return c;
}

void cohort_free(struct cohort *c)
{
	if (c) {
		cohort_deinit(c);
		xfree(c);
	}
}

void cohort_init(struct cohort *c, const struct vector *x)
{
	assert(c);
	assert(x);

	vector_init_copy(&c->traits, x);
	intset_init(&c->members);
}

void cohort_deinit(struct cohort *c)
{
	assert(c);
	intset_deinit(&c->members);
	vector_deinit(&c->traits);
}

const struct vector *cohort_traits(const struct cohort *c)
{
	assert(c);
	return &c->traits;
}

bool cohort_contains(const struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	assert(id <= INTPTR_MAX);
	return intset_contains(&c->members, id);
}

bool cohort_add(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	assert(id <= INTPTR_MAX);
	return intset_add(&c->members, id);
}

bool cohort_remove(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	assert(id <= INTPTR_MAX);
	return intset_remove(&c->members, id);
}

struct cohort_iter cohort_iter_make(const struct cohort *c)
{
	assert(c);

	struct cohort_iter it;
	it.members_it = intset_iter_make(&c->members);
	return it;
}

bool cohort_iter_advance(struct cohort_iter *it)
{
	assert(it);
	return intset_iter_advance(&it->members_it);
}

void cohort_iter_reset(struct cohort_iter *it)
{
	assert(it);
	intset_iter_reset(&it->members_it);
}
