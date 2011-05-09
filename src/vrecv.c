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
				     ssize_t index,
				     void *udata)
{
	assert(v);
	assert(v->dim >= 0);
	assert(e);
	assert(f);
	assert(f->design);	
	assert(index >= 0);
	assert(index <= design_dim(f->design) - v->dim);
	assert(udata);
	
	struct vrecv_frame *vf = udata;
	const struct vrecv_active *key = container_of(&e->dyad, const struct vrecv_active, dyad);
	struct vrecv_active *active;
	struct svector *dx = frame_dx(f, e->dyad.jrecv, e->dyad.isend);
	
	if (!dx)
		return false;
	
	if (e->type == DYAD_EVENT_INIT) {
		struct hashset_pos pos;
		
		if ((active = hashset_find(&vf->active, key, &pos)))
			goto move;
		
		if ((active = hashset_insert(&vf->active, &pos, key)))
			goto init;

		return false;
	} else { // e->type == DYAD_EVENT_MOVE
		active = hashset_lookup(&vf->active, key);
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

static bool vrecv_frame_init(struct vrecv_frame *vf)
{
	if (!hashset_init(&vf->active, vrecv_active_hash, vrecv_active_equals,
			  sizeof(struct vrecv_active)))
		goto fail_active;

	return true;
fail_active:
	return false;
}

static void *vrecv_frame_alloc(struct dyad_var *v, const struct frame *f)
{
	assert(v);
	assert(v->dim >= 0);
	assert(f);
	assert(f->design);	
	
	struct vrecv_frame *vf = malloc(sizeof(*vf));
	
	if (vf) {
		if (vrecv_frame_init(vf))
			return vf;
		free(vf);
	}
	return NULL;

}

static void vrecv_frame_deinit(struct vrecv_frame *vf)
{
	assert(vf);
	hashset_deinit(&vf->active);

}

static void vrecv_frame_free(struct dyad_var *v, const struct frame *f, void *udata)
{
	assert(v);
	assert(f);
	
	struct vrecv_frame *vf = udata;
	vrecv_frame_deinit(vf);
	free(vf);
}

static void vrecv_frame_clear(struct dyad_var *v, const struct frame *f, void *udata)
{
	assert(v);
	assert(f);
	
	struct vrecv_frame *vf = udata;
	hashset_clear(&vf->active);
}

bool vrecv_init(struct vrecv *v, const struct design *d)
{
	assert(v);
	assert(d);
	
	
	ssize_t n = vector_dim(design_intervals(d));
	v->dyad_var.dim = n + 1;
	v->dyad_var.dyad_event_mask = DYAD_EVENT_INIT | DYAD_EVENT_MOVE;
	v->dyad_var.handle_dyad_event = vrecv_handle_dyad_event;
	v->dyad_var.frame_alloc = vrecv_frame_alloc;
	v->dyad_var.frame_free = vrecv_frame_free;
	v->dyad_var.frame_clear = vrecv_frame_clear;
	return true;
	    
}

void vrecv_deinit(struct vrecv *v)
{
	assert(v);
}
