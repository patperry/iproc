#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vnrecv.h"

static void vnrecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
}

static void vnrecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static void vnrecv_handle_event(struct frame_var *fv,
				const struct frame_event *e, struct frame *f)
{
	assert(fv);
	assert(e);
	assert(f);
	assert(fv->design);
	assert(fv->design->index >= 0);
	assert(fv->design->index <= design_dim(f->design) - fv->design->dim);
	assert(e->type & (DYAD_EVENT_INIT | DYAD_EVENT_MOVE));

	const struct dyad_event_meta *meta = &e->meta.dyad;
	struct frame_event dx;

	dx.type = DYAD_VAR_EVENT;
	dx.time = e->time;
	dx.id = -1;
	dx.meta.dyad_var.item.isend = meta->msg_dyad.jrecv;
	dx.meta.dyad_var.item.jrecv = meta->msg_dyad.isend;

	if (e->type == DYAD_EVENT_MOVE) {
		dx.meta.dyad_var.index = fv->design->index + meta->intvl - 1;
		dx.meta.dyad_var.delta = -1.0;
		frame_events_add(f, &dx);
	}

	dx.meta.dyad_var.index = fv->design->index + meta->intvl;
	dx.meta.dyad_var.delta = 1.0;
	frame_events_add(f, &dx);
}

static struct var_type VAR_TYPE_NRECV_REP = {
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	vnrecv_init,
	vnrecv_deinit,
	NULL,			// frame_init,
	NULL,			// frame_deinit,
	NULL,			// frame_clear,
	vnrecv_handle_event
};

const struct var_type *VAR_TYPE_NRECV = &VAR_TYPE_NRECV_REP;
