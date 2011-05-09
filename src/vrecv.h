#ifndef _VRECV
#define _VRECV

#include "design.h"
#include "hashset.h"
#include "intmap.h"

struct vrecv_active {
	struct dyad dyad;
	ssize_t id;
	ssize_t intvl;
};

struct vrecv_frame {
	struct hashset active;
};

struct vrecv {
	struct dyad_var dyad_var;
};

bool vrecv_init(struct vrecv *v, const struct design *d);
void vrecv_deinit(struct vrecv *v);

#endif /* _VRECV */
