#include "port.h"
#include <assert.h>
#include "strata.h"

void strata_init(struct strata *s, ssize_t dim)
{
	assert(s);
	assert(dim >= 0);

	s->dim = dim;
	array_init(&s->levels, sizeof(struct vector));
	intmap_init(&s->level_hashes, sizeof(struct intset),
		    alignof(struct intset));
}

void strata_deinit(struct strata *s)
{
	assert(s);

	struct intmap_iter it;
	INTMAP_FOREACH(it, &s->level_hashes) {
		struct intset *set = INTMAP_VAL(it);
		intset_deinit(set);
	}
	intmap_deinit(&s->level_hashes);

	struct vector *level;
	ARRAY_FOREACH(level, &s->levels) {
		vector_deinit(level);
	}
	array_deinit(&s->levels);
}

static ssize_t lookup_add(struct strata *s, const struct vector *level,
			  bool add)
{
	assert(s);
	assert(level);
	assert(vector_dim(level) == strata_dim(s));

	struct vector *stratum;
	ssize_t index;

	int32_t hash = vector_hash(level);
	struct intmap_pos pos;
	struct intset *set = intmap_find(&s->level_hashes, hash, &pos);

	if (!set) {
		set = intmap_insert(&s->level_hashes, &pos, NULL);
		intset_init(set);
	}

	struct intset_iter it;
	INTSET_FOREACH(it, set) {
		index = INTSET_KEY(it);
		stratum = array_item(&s->levels, index);
		if (vector_equals(stratum, level)) {
			goto found;
		}
	}
	index = -1;		/* level not found */

	/* level not found; create a new stratum */
	if (add) {
		index = array_count(&s->levels);
		stratum = array_add(&s->levels, NULL);
		vector_init_copy(stratum, level);
		intset_add(set, index);
	}
found:
	return index;
}

ssize_t strata_add(struct strata *s, const struct vector *level)
{
	assert(s);
	assert(level);
	assert(vector_dim(level) == strata_dim(s));

	return lookup_add(s, level, true);
}

ssize_t strata_lookup(const struct strata *s, const struct vector *level)
{
	assert(s);
	assert(level);
	assert(vector_dim(level) == strata_dim(s));

	return lookup_add((struct strata *)s, level, false);
}
