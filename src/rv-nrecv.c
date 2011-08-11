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

static void nrecv_init(struct design_var *dv, const struct design *d, void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

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

static void handle_clear(void *udata, struct frame *f)
{
	struct frame_var *fv = udata;
	struct nrecv_udata *fv_udata = fv->udata;
	
	assert(fv);
	assert(fv->udata);

	struct nrecv_active *active;

	ARRAY_FOREACH(active, &fv_udata->active) {
		intset_clear(&active->jrecv);
	}
}

static void handle_message_add(void *udata, struct frame *f, const struct message *msg)
{
	struct frame_var *fv = udata;
	struct nrecv_udata *fv_udata = fv->udata;	
	
	assert(fv);
	assert(f);
	assert(msg);
	assert(fv->design);
	assert(fv->design->dyn_index >= 0);
	assert(fv->design->dyn_index + fv->design->dim
	       <= design_recv_dyn_dim(f->design));

	
	ssize_t jrecv = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t isend = msg->to[ito];

		struct nrecv_active *active = array_item(&fv_udata->active, isend);
		if (intset_add(&active->jrecv, jrecv)) {
			frame_recv_update(f, isend, jrecv, dyn_index, 1.0);
		}

		
		ssize_t i, n = fv->design->dim;
		for (i = 1; i < n; i++) {
			frame_recv_update(f, isend, jrecv, dyn_index + i, 1.0);
		}
	}
}

static void handle_message_advance(void *udata, struct frame *f, const struct message *msg, ssize_t intvl)
{
	struct frame_var *fv = udata;
	
	assert(fv);
	assert(f);
	assert(msg);
	assert(fv->design);
	assert(fv->design->dyn_index >= 0);
	assert(fv->design->dyn_index + fv->design->dim
	       <= design_recv_dyn_dim(f->design));
		
	ssize_t jrecv = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t isend = msg->to[ito];
		
		frame_recv_update(f, isend, jrecv, dyn_index + intvl, -1.0);
	}
}


static struct var_type RECV_VAR_NRECV_REP = {
	VAR_RECV_VAR,
	nrecv_init,
	nrecv_deinit,
	nrecv_frame_init,
	nrecv_frame_deinit,
	{
		handle_message_add,
		handle_message_advance,
		NULL, // recv_update
		NULL, // send_update
		handle_clear
	}
};

const struct var_type *RECV_VAR_NRECV = &RECV_VAR_NRECV_REP;
