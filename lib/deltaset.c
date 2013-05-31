#include "port.h"
#include "xalloc.h"
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "deltaset.h"

#ifdef DEBUG_DELTASET

static void deltaset_check(const struct deltaset *ds)
{
	const struct delta *prev, *cur, *next;

	size_t i, n = deltaset_count(ds);
	for (i = 0; i < n; i++) {
		cur = &ds->delta[i];
		assert(cur->i == i);
		prev = cur->prev;
		next = cur->next;

		if (prev) {
			assert(prev->next == cur);
			assert(cur->t >= prev->t);
		}
		if (next) {
			assert(next->prev == cur);
			assert(cur->t <= next->t);
		} else {
			assert(ds->head == cur);
		}
	}
}

#define CHECK(ds) deltaset_check(ds)

#else

#define CHECK(ds) (void)(ds)

#endif




void deltaset_init(struct deltaset *ds, size_t n)
{
	ds->n = n;
	ds->delta = xmalloc(n * sizeof(*ds->delta));
	deltaset_clear(ds);
	CHECK(ds);
}


void deltaset_deinit(struct deltaset *ds)
{
	free(ds->delta);
}


void deltaset_clear(struct deltaset *ds)
{
	size_t i, n = deltaset_count(ds);
	struct delta *prev, *cur, *next;

	ds->head = n ? &ds->delta[0] : NULL;

	prev = ds->head;
	cur = NULL;
	for (i = 0; i < n; i++) {
		next = cur;
		cur = prev;
		prev = i + 1 < n ? &ds->delta[i+1] : NULL;

		cur->i = i;
		cur->t = -INFINITY;
		cur->prev = prev;
		cur->next = next;
	}
	CHECK(ds);
}


void deltaset_update(struct deltaset *ds, size_t i, double t)
{
	// printf("update(%p, %zd, %f)\n", ds, i, t);

	assert(i < deltaset_count(ds));
	assert(t >= deltaset_tlast(ds, i));
	assert(!isnan(t));

	struct delta *cur = &ds->delta[i];

	assert(cur->i == i);

	/* remove delta[i] from the update list */
	struct delta *prev = cur->prev;
	struct delta *next = cur->next;

	assert(!next || next->prev == cur);
	assert(!prev || prev->next == cur);
	if (next) {
		next->prev = prev;
	} else {
		ds->head = prev;
	}
	if (prev)
		prev->next = next;


	/* find appropriate place for delta[i] */
	if (!ds->head || ds->head->t <= t) {	/* delta[i] belongs at new head */
		cur->prev = ds->head;
		cur->next = NULL;
		ds->head->next = cur;
		ds->head = cur;
	} else {
		next = ds->head;
		prev = ds->head->prev;
		
		while (prev && prev->t > t) {
			assert(t < next->t);
			next = prev;
			prev = next->prev;
		}

		/* delta[i] belongs between next and prev */
		assert(next);
		next->prev = cur;
		cur->next = next;

		if (prev) {
			prev->next = cur;
		}
		cur->prev = prev;
	}

	/* update time */
	cur->t = t;

	CHECK(ds);
}
