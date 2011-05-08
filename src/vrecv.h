#ifndef _VRECV
#define _VRECV

#include "design.h"
#include "hashset.h"

struct vrecv_active {
	struct dyad dyad;
	ssize_t id;
	ssize_t intvl;
};

struct vrecv {
	struct dyad_var dyad_var;
	struct hashset active;
};

bool vrecv_init(struct vrecv *v, const struct design *d);
void vrecv_deinit(struct vrecv *v);

#endif /* _VRECV */
