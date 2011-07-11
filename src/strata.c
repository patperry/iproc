#include "port.h"
#include <assert.h>
#include "strata.h"

static void stratum_init(struct stratum *s, const struct vector *traits)
{
	assert(s);
	assert(traits);
	vector_init_copy(&s->traits, traits);
	intset_init(&s->ids);
}


static void stratum_deinit(struct stratum *s)
{
	assert(s);
	intset_deinit(&s->ids);
	vector_deinit(&s->traits);
}


static bool stratum_add(struct stratum *s, intptr_t id)
{
	assert(s);
	return intset_add(&s->ids, id);
}


void strata_init(struct strata *s, ssize_t dim)
{
	assert(s);
	assert(dim >= 0);
	
	s->dim = dim;
	array_init(&s->items, sizeof(struct stratum));
	intmap_init(&s->trait_hashes, sizeof(struct intset),
		    alignof(struct intset));
}


void strata_deinit(struct strata *s)
{
	assert(s);
	
	struct intmap_iter it;
	INTMAP_FOREACH(it, &s->trait_hashes) {
		struct intset *set = INTMAP_VAL(it);
		intset_deinit(set);
	}
	intmap_deinit(&s->trait_hashes);
	
	struct stratum *item;
	ARRAY_FOREACH(item, &s->items) {
		stratum_deinit(item);
	}
	array_deinit(&s->items);
}


ssize_t strata_add(struct strata *s, intptr_t id, const struct vector *traits)
{
	assert(s);
	assert(vector_dim(traits) == strata_dim(s));

	struct stratum *item;
	ssize_t index = -1;

	int32_t hash = vector_hash(traits);
	struct intmap_pos pos;
	struct intset *set = intmap_find(&s->trait_hashes, hash, &pos);
	
	if (!set) {
		set = intmap_insert(&s->trait_hashes, &pos, NULL);
		intset_init(set);
	}
	
	struct intset_iter it;
	INTSET_FOREACH(it, set) {
		index = INTSET_KEY(it);
		item = array_item(&s->items, index);
		const struct vector *x = stratum_traits(item);
		if (vector_equals(x, traits)) {
			goto found;
		}
	}
	
	/* cohort not found; create a new one */
	index = array_count(&s->items);
	item = array_add(&s->items, NULL);
	stratum_init(item, traits);
	intset_add(set, index);
found:
	stratum_add(item, id);
	return index;
}
