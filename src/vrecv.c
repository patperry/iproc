#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "design.h"
#include "util.h"
#include "vrecv.h"
#include "frame.h"

uint32_t vrecv_active_hash (const void *x)
{
	const struct vrecv_active *active = x;
	return memory_hash(&active->dyad, sizeof(active->dyad));
}

bool vrecv_active_equals (const void *x, const void *y)
{
	const struct vrecv_active *a = x, *b = y;
	return a->dyad.isend == b->dyad.isend && a->dyad.jrecv == b->dyad.jrecv;
}

static bool vrecv_handle_dyad_event (struct dyad_var *v,
				     const struct dyad_event *e,
				     struct frame *f,
				     ssize_t index)
{
	assert(v);
	assert(v->dim >= 0);
	assert(e);
	assert(f);
	assert(f->design);	
	assert(index >= 0);
	assert(index <= design_dim(f->design) - v->dim);
	

	struct vrecv *vrecv = container_of(v, struct vrecv, dyad_var);
	const struct vrecv_active *key = container_of(&e->dyad, const struct vrecv_active, dyad);
	struct vrecv_active *active;
	struct svector *dx = frame_dx(f, e->dyad.jrecv, e->dyad.isend);
	
	if (!dx)
		return false;
	
	if (e->type == DYAD_EVENT_INIT) {
		struct hashset_pos pos;
		
		if ((active = hashset_find(&vrecv->active, key, &pos)))
			goto move;
		
		if ((active = hashset_insert(&vrecv->active, &pos, key)))
			goto init;

		return false;
	} else { // e->type == DYAD_EVENT_MOVE
		active = hashset_lookup(&vrecv->active, key);
		assert(active);

		if (active->id == e->id) {
			assert(active->intvl == e->intvl - 1);
			goto move;
		}
		
		return true;
	}

move:		
	svector_clear_at(dx, index + active->intvl);
init:
	if (!svector_set(dx, index + e->intvl, 1.0))
		return false;
	active->id = e->id;
	active->intvl = e->intvl;
	return true;
}


bool vrecv_init(struct vrecv *v, const struct design *d)
{
	assert(v);
	assert(d);
	
	if (!hashset_init(&v->active, vrecv_active_hash, vrecv_active_equals,
			  sizeof(struct vrecv_active)))
	    goto fail_active;
	
	ssize_t n = vector_dim(design_intervals(d));
	v->dyad_var.dim = n + 1;
	v->dyad_var.dyad_event_mask = DYAD_EVENT_INIT | DYAD_EVENT_MOVE;
	v->dyad_var.handle_dyad_event = vrecv_handle_dyad_event;
	return true;
	    
fail_active:
	return false;
}

void vrecv_deinit(struct vrecv *v)
{
	assert(v);
	hashset_deinit(&v->active);
}
