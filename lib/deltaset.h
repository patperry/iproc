#ifndef DELTASET_H
#define DELTASET_H

#include <assert.h>
#include <math.h>
#include <stddef.h>


struct delta {
	size_t i;
	double t;
	struct delta *prev, *next;
};

#define DELTA_NEXT(d) ((d)->next)
#define DELTA_PREV(d) ((d)->prev)
#define DELTA_ITEM(d) ((d)->i)
#define DELTA_TIME(d) ((d)->t)


#define DELTASET_FOREACH(d,ds) \
	for (d = deltaset_head(ds); d != NULL; d = DELTA_PREV(d))


struct deltaset {
	size_t n;
	struct delta *delta;
	struct delta *head;
};



void deltaset_init(struct deltaset *ds, size_t n);
void deltaset_deinit(struct deltaset *ds);
void deltaset_clear(struct deltaset *ds);

static inline size_t deltaset_count(const struct deltaset *ds);
static inline double deltaset_tlast(const struct deltaset *ds, size_t i);
static inline double deltaset_thead(const struct deltaset *ds);
static inline struct delta *deltaset_head(const struct deltaset *ds);


void deltaset_update(struct deltaset *ds, size_t i, double t);


/* inline function definitions */
size_t deltaset_count(const struct deltaset *ds)
{
	return ds->n;
}

double deltaset_tlast(const struct deltaset *ds, size_t i)
{
	assert(i < deltaset_count(ds));
	return ds->delta[i].t;
}

double deltaset_thead(const struct deltaset *ds)
{
	return ds->head ? ds->head->t : -INFINITY;
}

struct delta *deltaset_head(const struct deltaset *ds){
	return ds->head;
}



#endif /* DELTASET_H */