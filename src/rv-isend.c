#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void isend_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

	dv->dim = 1;
}

static void isend_message_add(void *udata, struct frame *f,
			       const struct message *msg)
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
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		if (vector_item(dx, dyn_index) == 0.0) {
			frame_recv_update(f, isend, jrecv, dyn_index, 1.0);
		}
	}
}

static struct var_type RECV_VAR_ISEND_REP = {
	VAR_RECV_VAR,
	isend_init,
	NULL, // deinit
	NULL, // frame_init
	NULL, // frame_deinit
	{
		isend_message_add,
		NULL,			// message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL			// clear
	}
};

const struct var_type *RECV_VAR_ISEND = &RECV_VAR_ISEND_REP;
