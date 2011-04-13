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

bool cohort_add(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	return intset_add(&c->members, id);
}

void cohort_remove(struct cohort *c, ssize_t id)
{
	assert(c);
	assert(id >= 0);
	intset_remove(&c->members, id);
}
