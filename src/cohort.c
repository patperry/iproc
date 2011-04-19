#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "cohort.h"

struct cohort *cohort_new(const struct vector *x)
{
	assert(x);

	struct cohort *c = malloc(sizeof(*c));

	if (c) {
		if (cohort_init(c, x)) {
			return c;
		}
		free(c);
	}
	return NULL;
}

void cohort_free(struct cohort *c)
{
	if (c) {
		cohort_deinit(c);
		free(c);
	}
}

bool cohort_init(struct cohort *c, const struct vector *x)
{
	assert(c);
	assert(x);

	if (vector_init_copy(&c->traits, x)) {
		if (intset_init(&c->members)) {
			return c;
		}
		vector_deinit(&c->traits);
	}
	return NULL;
}

void cohort_deinit(struct cohort *c)
{
	assert(c);
	intset_deinit(&c->members);
	vector_deinit(&c->traits);
}

ssize_t cohort_dim(const struct cohort *c)
{
	assert(c);
	return vector_size(cohort_traits(c));
}

const struct vector *cohort_traits(const struct cohort *c)
{
	assert(c);
	return &c->traits;
}

bool cohort_empty(const struct cohort * c)
{
	assert(c);
	return intset_empty(&c->members);
}

ssize_t cohort_size(const struct cohort * c)
{
	assert(c);
	return intset_size(&c->members);
}

bool cohort_contains(const struct cohort * c, ssize_t id)
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

void cohort_remove(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	assert(id <= INTPTR_MAX);
	intset_remove(&c->members, id);
}

void cohort_iter_init(const struct cohort *c, struct cohort_iter *it)
{
	assert(c);
	assert(it);
	intset_iter_init(&c->members, &it->members_it);
}

void cohort_iter_deinit(const struct cohort *c, struct cohort_iter *it)
{
	assert(c);
	assert(it);
	intset_iter_deinit(&c->members, &it->members_it);
}

bool cohort_iter_advance(const struct cohort *c, struct cohort_iter *it)
{
	assert(c);
	assert(it);
	return intset_iter_advance(&c->members, &it->members_it);
}

ssize_t cohort_iter_current(const struct cohort *c,
			    const struct cohort_iter *it)
{
	assert(c);
	assert(it);
	return (ssize_t)intset_iter_current(&c->members, &it->members_it);
}
