#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void nsend_init(struct design_var *dv, const struct design *d, void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
}

static void nsend_deinit(struct design_var *dv)
{
	assert(dv);
}

static void nsend_frame_init(struct frame_var *fv, struct frame *f)
{
	assert(fv);
	assert(f);
}

static void nsend_frame_deinit(struct frame_var *fv)
{
	assert(fv);
}

static void handle_message_add(void *udata, struct frame *f, const struct message *msg)
{
	struct frame_var *fv = udata;
	
	assert(fv);
	assert(f);
	assert(msg);
	assert(fv->design);
	assert(fv->design->dyn_index >= 0);
	assert(fv->design->dyn_index + fv->design->dim
	       <= design_recv_dyn_dim(f->design));

	
	ssize_t isend = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t jrecv = msg->to[ito];

		ssize_t i, n = fv->design->dim;
		for (i = 0; i < n; i++) {
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
		
	ssize_t isend = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t jrecv = msg->to[ito];
		
		frame_recv_update(f, isend, jrecv, dyn_index + intvl - 1, -1.0);
	}
}


static struct var_type RECV_VAR_NSEND_REP = {
	VAR_RECV_VAR,
	nsend_init,
	nsend_deinit,
	nsend_frame_init,
	nsend_frame_deinit,
	{
		handle_message_add,
		handle_message_advance,
		NULL, // recv_update
		NULL, // send_update
		NULL, // clear
	}
};

const struct var_type *RECV_VAR_NSEND = &RECV_VAR_NSEND_REP;
