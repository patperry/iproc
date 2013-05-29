#include "port.h"
#include "xalloc.h"
#include <assert.h>
#include <stdlib.h>
#include "deltaset.h"

void deltaset_init(struct deltaset *ds, size_t n)
{
	ds->n = n;
	ds->delta = xmalloc(n * sizeof(*ds->delta));
	deltaset_clear(ds);
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
}


void deltaset_update(struct deltaset *ds, size_t i, double t)
{
	assert(i < deltaset_count(ds));
	assert(t >= ds->delta[i].t);
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
		prev = ds->head->next;
		
		while (prev && prev->t < t) {
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
}
