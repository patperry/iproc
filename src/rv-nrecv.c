#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

struct nrecv_active {
	struct intset jrecv;
};

struct nrecv_udata {
	struct array active;
};


static void nrecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 2;
}

static void nrecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static void nrecv_frame_init(struct frame_var *fv, struct frame *f)
{
	assert(fv);
	assert(f);
	
	const struct design *d = frame_design(f);
	struct nrecv_udata *udata = xcalloc(1, sizeof(*udata));
	array_init(&udata->active, sizeof(struct nrecv_active));
	
	ssize_t i, n = design_recv_count(d);
	for (i = 0; i < n; i++) {
		struct nrecv_active *active = array_add(&udata->active, NULL);
		intset_init(&active->jrecv);
	}
	
	fv->udata = udata;
}

static void nrecv_frame_deinit(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);
	
	struct nrecv_udata *udata = fv->udata;
	struct nrecv_active *active;

	ARRAY_FOREACH(active, &udata->active) {
		intset_deinit(&active->jrecv);
	}
	
	array_deinit(&udata->active);
	xfree(udata);
}

static void nrecv_frame_clear(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);
	
	struct nrecv_udata *udata = fv->udata;
	struct nrecv_active *active;
	
	ARRAY_FOREACH(active, &udata->active) {
		intset_clear(&active->jrecv);
	}
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
	
	if (e->type == DYAD_EVENT_INIT) {
		struct nrecv_udata *udata = fv->udata;
		struct nrecv_active *active = array_item(&udata->active, meta->msg_dyad.jrecv);
		if (intset_add(&active->jrecv, meta->msg_dyad.isend)) {
			dx.type = RECV_VAR_EVENT;
			dx.time = e->time;
			dx.id = -1;
			dx.meta.recv_var.item.isend = meta->msg_dyad.jrecv;
			dx.meta.recv_var.item.jrecv = meta->msg_dyad.isend;
			dx.meta.recv_var.index = fv->design->index;
			dx.meta.recv_var.delta = 1.0;
			frame_events_add(f, &dx);
		}
	}

	dx.type = RECV_VAR_EVENT;
	dx.time = e->time;
	dx.id = -1;
	dx.meta.recv_var.item.isend = meta->msg_dyad.jrecv;
	dx.meta.recv_var.item.jrecv = meta->msg_dyad.isend;

	if (e->type == DYAD_EVENT_MOVE) {
		assert(meta->intvl > 0);
		dx.meta.recv_var.index = fv->design->index + meta->intvl;
		dx.meta.recv_var.delta = -1.0;
		frame_events_add(f, &dx);
	}

	dx.meta.recv_var.index = fv->design->index + 1 + meta->intvl;
	dx.meta.recv_var.delta = 1.0;
	frame_events_add(f, &dx);
}

static struct var_type RECV_VAR_NRECV_REP = {
	VAR_RECV_VAR,
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	nrecv_init,
	nrecv_deinit,
	nrecv_frame_init,
	nrecv_frame_deinit,
	nrecv_frame_clear,
	nrecv_handle_event
};

const struct var_type *RECV_VAR_NRECV = &RECV_VAR_NRECV_REP;
