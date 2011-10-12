#include "port.h"
#include <stdint.h>		// uint64_t
#include <stdlib.h>		// free
#include <string.h>		// memcpy
#include "hash.h"		// double_hash, hash_combine
#include "util.h"		// container_of
#include "xalloc.h"		// xmalloc, xrealloc
#include "strata.h"

struct level {
	double *data;
	size_t index;
};

static size_t level_hash(const struct hashset *set, const void *x)
{
	const double *level = *(const double **)x;
	size_t i, n = container_of(set, struct strata, levels_set)->dim;
	size_t seed = 0;

	for (i = 0; i < n; i++) {
		seed = hash_combine(seed, double_hash(level[i]));
	}

	return seed;
}

static int level_compar(const struct hashset *set, const void *x, const void *y)
{
	const uint64_t *u = *(const uint64_t **)x;
	const uint64_t *v = *(const uint64_t **)y;
	size_t i, n = container_of(set, struct strata, levels_set)->dim;

	for (i = 0; i < n; i++) {
		if (u[i] != v[i])
			return 1;
	}

	return 0;
}

void strata_init(struct strata *s, size_t dim)
{
	s->dim = dim;
	s->nlevel = 0;
	s->nlevel_max = 0;
	s->levels = NULL;
	hashset_init(&s->levels_set, sizeof(struct level), level_hash,
		     level_compar);
}

void strata_deinit(struct strata *s)
{
	strata_clear(s);
	hashset_deinit(&s->levels_set);
	free(s->levels);
}

void strata_clear(struct strata *s)
{
	hashset_clear(&s->levels_set);

	double **l = s->levels;
	size_t i, n = s->nlevel;

	for (i = 0; i < n; i++) {
		free(l[i]);
	}

	s->nlevel = 0;
}

size_t strata_add(struct strata *s, const double *x)
{
	struct hashset_pos pos;
	struct level *l;

	if ((l = hashset_find(&s->levels_set, &x, &pos))) {
		return l->index;
	}

	/* expand the levels array if necessary */
	if (s->nlevel >= s->nlevel_max) {
		if (!s->nlevel_max) {
			s->nlevel_max = 5;
		} else {
			s->nlevel_max *= 2;
		}
		s->levels = xrealloc(s->levels,
				     s->nlevel_max * sizeof(s->levels[0]));
	}

	/* duplicate x */
	size_t n = s->dim;
	double *x1 = xmalloc(n * sizeof(x[0]));
	memcpy(x1, x, n * sizeof(x[0]));

	/* insert the duplicate into the levels array */
	size_t index = s->nlevel;
	s->levels[index] = x1;
	s->nlevel++;

	/* insert the new level into the levels set */
	struct level l1 = { x1, index };
	hashset_insert(&s->levels_set, &pos, &l1);

	return index;
}

ptrdiff_t strata_find(const struct strata *s, const double *level)
{
	const struct level *l = hashset_item(&s->levels_set, &level);
	if (l) {
		return l->index;
	} else {
		return -1;
	}
}
