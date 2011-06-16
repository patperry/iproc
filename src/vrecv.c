#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "design.h"
#include "frame.h"
#include "hashset.h"
#include "intmap.h"
#include "util.h"
#include "vrecv.h"

struct vrecv_active {
	struct dyad dyad;
	ssize_t id;
	ssize_t intvl;
};

struct vrecv_udata {
	struct hashset active;
};

static uint32_t vrecv_active_hash(const void *x)
{
	const struct vrecv_active *active = x;
	return memory_hash(&active->dyad, sizeof(active->dyad));
}

static bool vrecv_active_equals(const void *x, const void *y)
{
	const struct vrecv_active *a = x, *b = y;
	return a->dyad.isend == b->dyad.isend && a->dyad.jrecv == b->dyad.jrecv;
}

static void vrecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
}

static void vrecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static void vrecv_frame_init(struct frame_var *fv, struct frame *f)
{
	assert(fv);
	assert(f);

	struct vrecv_udata *udata = xcalloc(1, sizeof(*udata));
	hashset_init(&udata->active, vrecv_active_hash, vrecv_active_equals,
		     sizeof(struct vrecv_active));
	fv->udata = udata;
}

static void vrecv_frame_deinit(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);

	struct vrecv_udata *udata = fv->udata;
	hashset_deinit(&udata->active);
	free(udata);
}

static void vrecv_frame_clear(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);

	struct vrecv_udata *udata = fv->udata;
	hashset_clear(&udata->active);
}

static void vrecv_handle_event(struct frame_var *fv,
			       const struct frame_event *e, struct frame *f)
{
	assert(fv);
	assert(e);
	assert(f);
	assert(fv->design);
	assert(fv->design->index >= 0);
	assert(fv->design->index <=
	       design_recv_dim(f->design) - fv->design->dim);
	assert(fv->udata);
	assert(e->type & (DYAD_EVENT_INIT | DYAD_EVENT_MOVE));

	const struct dyad_event_meta *meta = &e->meta.dyad;

	ssize_t index = fv->design->index;
	struct vrecv_udata *udata = fv->udata;
	const struct vrecv_active *key =
	    container_of(&meta->msg_dyad, const struct vrecv_active, dyad);
	struct vrecv_active *active;

	struct frame_event dx;
	dx.type = RECV_VAR_EVENT;
	dx.time = e->time;
	dx.id = -1;
	dx.meta.recv_var.item.isend = meta->msg_dyad.jrecv;
	dx.meta.recv_var.item.jrecv = meta->msg_dyad.isend;

	if (e->type == DYAD_EVENT_INIT) {
		struct hashset_pos pos;

		if ((active = hashset_find(&udata->active, key, &pos)))
			goto move;

		active = hashset_insert(&udata->active, &pos, key);
		goto init;
	} else {		// e->type == DYAD_EVENT_MOVE
		active = hashset_item(&udata->active, key);
		assert(active);

		if (active->id == e->id) {
			assert(active->intvl == meta->intvl - 1);
			goto move;
		}
		goto out;
	}
move:
	dx.meta.recv_var.index = index + active->intvl;
	dx.meta.recv_var.delta = -1.0;
	frame_events_add(f, &dx);
init:
	dx.meta.recv_var.index = index + meta->intvl;
	dx.meta.recv_var.delta = +1.0;
	frame_events_add(f, &dx);

	active->id = e->id;
	active->intvl = meta->intvl;
out:
	return;
}

static struct var_type VAR_TYPE_RECV_REP = {
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	vrecv_init,
	vrecv_deinit,
	vrecv_frame_init,
	vrecv_frame_deinit,
	vrecv_frame_clear,
	vrecv_handle_event
};

const struct var_type *VAR_TYPE_RECV = &VAR_TYPE_RECV_REP;
