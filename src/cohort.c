#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "compare.h"
#include "cohort.h"

DEFINE_EQUALS_FN(ssize_equals, ssize_t)


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
		if (darray_init(&c->members, sizeof(ssize_t))) {
			return c;
		}
		vector_deinit(&c->traits);
	}
	return NULL;
}

void cohort_deinit(struct cohort *c)
{
	assert(c);
	darray_deinit(&c->members);
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

bool cohort_empty(const struct cohort *c)
{
	assert(c);
	return darray_empty(&c->members);
}

ssize_t cohort_size(const struct cohort *c)
{
	assert(c);
	return darray_size(&c->members);
}

bool cohort_contains(const struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	return darray_contains(&c->members, &id, ssize_equals);
}

bool cohort_add(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	if (!cohort_contains(c, id)) {
		return darray_push_back(&c->members, &id);
	}
	return true;
}

void cohort_remove(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);

	ssize_t i = darray_find_last_index(&c->members, &id, ssize_equals);
	if (i >= 0)
		darray_erase(&c->members, i);
}

void cohort_iter_init(const struct cohort *c, struct cohort_iter *it)
{
	assert(c);
	assert(it);
	it->index = -1;
}

void cohort_iter_deinit(const struct cohort *c, struct cohort_iter *it)
{
	assert(c);
	assert(it);
}

bool cohort_iter_advance(const struct cohort *c, struct cohort_iter *it)
{
	assert(c);
	assert(it);
	it->index++;
	return it->index < darray_size(&c->members);
}

ssize_t cohort_iter_current(const struct cohort *c,
			    const struct cohort_iter *it)
{
	assert(c);
	assert(it);
	return *(ssize_t *)darray_at(&c->members, it->index);
}
