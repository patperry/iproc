#ifndef VERSION_H
#define VERSION_H

#include <stddef.h>


struct version {
	size_t counter;
	struct version_watch **w;
	size_t nw, nwmax;
};

struct version_watch {
	size_t counter;
};


void version_init(struct version *v);
void version_deinit(struct version *v);

void version_update(struct version *v);


static inline int version_changed(const struct version *v,
				  const struct version_watch *w)
{
	return w->counter < v->counter;
}


void version_watch_init(struct version_watch *w, struct version *v);
void version_watch_deinit(struct version_watch *w, struct version *v);

static inline void version_watch_set(struct version_watch *w,
				     const struct version *v)
{
	w->counter = v->counter;
}



#endif /* VERSION_H */
