#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void nrecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
}

static void nrecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static void nrecv_handle_event(struct frame_var *fv,
			       const struct frame_event *e, struct frame *f)
{
	assert(fv);
	assert(e);
	assert(f);
	assert(fv->design);
	assert(fv->design->index >= 0);
	assert(fv->design->index <=
	       design_recv_dim(f->design) - fv->design->dim);
	assert(e->type & (DYAD_EVENT_INIT | DYAD_EVENT_MOVE));

	const struct dyad_event_meta *meta = &e->meta.dyad;
	struct frame_event dx;

	dx.type = RECV_VAR_EVENT;
	dx.time = e->time;
	dx.id = -1;
	dx.meta.recv_var.item.isend = meta->msg_dyad.jrecv;
	dx.meta.recv_var.item.jrecv = meta->msg_dyad.isend;

	if (e->type == DYAD_EVENT_MOVE) {
		dx.meta.recv_var.index = fv->design->index + meta->intvl - 1;
		dx.meta.recv_var.delta = -1.0;
		frame_events_add(f, &dx);
	}

	dx.meta.recv_var.index = fv->design->index + meta->intvl;
	dx.meta.recv_var.delta = 1.0;
	frame_events_add(f, &dx);
}

static struct var_type RECV_VAR_NRECV_REP = {
	VAR_RECV_VAR,
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	nrecv_init,
	nrecv_deinit,
	NULL,			// frame_init,
	NULL,			// frame_deinit,
	NULL,			// frame_clear,
	nrecv_handle_event
};

const struct var_type *RECV_VAR_NRECV = &RECV_VAR_NRECV_REP;
