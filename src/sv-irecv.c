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
#include "vars.h"

struct irecv_active {
	ssize_t isend;
	ssize_t id;
	ssize_t intvl;
};

struct irecv_udata {
	struct hashset active;
};

static uint32_t irecv_active_hash(const void *x)
{
	const struct irecv_active *active = x;
	return memory_hash(&active->isend, sizeof(active->isend));
}

static bool irecv_active_equals(const void *x, const void *y)
{
	const struct irecv_active *a = x, *b = y;
	return a->isend == b->isend;
}

static void irecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
}

static void irecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static void irecv_frame_init(struct frame_var *fv, struct frame *f)
{
	assert(fv);
	assert(f);

	struct irecv_udata *udata = xcalloc(1, sizeof(*udata));
	hashset_init(&udata->active, irecv_active_hash, irecv_active_equals,
		     sizeof(struct irecv_active));
	fv->udata = udata;
}

static void irecv_frame_deinit(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);

	struct irecv_udata *udata = fv->udata;
	hashset_deinit(&udata->active);
	free(udata);
}

static void irecv_frame_clear(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);

	struct irecv_udata *udata = fv->udata;
	hashset_clear(&udata->active);
}

static void irecv_handle_event(struct frame_var *fv,
			       const struct frame_event *e, struct frame *f)
{
	assert(fv);
	assert(e);
	assert(f);
	assert(fv->design);
	assert(fv->design->index >= 0);
	assert(fv->design->index <=
	       design_send_dim(f->design) - fv->design->dim);
	assert(fv->udata);
	assert(e->type & (DYAD_EVENT_INIT | DYAD_EVENT_MOVE));

	const struct dyad_event_meta *meta = &e->meta.dyad;

	ssize_t index = fv->design->index;
	struct irecv_udata *udata = fv->udata;
	const struct irecv_active *key =
	    container_of(&meta->msg_dyad.jrecv, const struct irecv_active, isend);
	struct irecv_active *active;

	struct frame_event dx;
	dx.type = SEND_VAR_EVENT;
	dx.time = e->time;
	dx.id = -1;
	dx.meta.send_var.item = meta->msg_dyad.jrecv;

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
	dx.meta.send_var.index = index + active->intvl;
	dx.meta.send_var.delta = -1.0;
	frame_events_add(f, &dx);
init:
	dx.meta.send_var.index = index + meta->intvl;
	dx.meta.send_var.delta = +1.0;
	frame_events_add(f, &dx);

	active->id = e->id;
	active->intvl = meta->intvl;
out:
	return;
}

static struct var_type SEND_VAR_IRECV_REP = {
	VAR_SEND_VAR,
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	irecv_init,
	irecv_deinit,
	irecv_frame_init,
	irecv_frame_deinit,
	irecv_frame_clear,
	irecv_handle_event
};

const struct var_type *SEND_VAR_IRECV = &SEND_VAR_IRECV_REP;
