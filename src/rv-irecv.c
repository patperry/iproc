#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static char *irecv_names[] = { "IRecv" };

static void irecv_init(struct design_var *v, const struct design *d,
		       void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = irecv_names;
}

static void irecv_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t jrecv = msg->from;
	size_t dyn_index = v->dyn_index;
	double dx_data[1] = { 1.0 };
	size_t dx_index[1] = { dyn_index };
	size_t dx_nnz = 1;

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		if (vector_item(dx, dyn_index) == 0.0) {
			frame_recv_update(f, isend, jrecv, dx_data, dx_index,
					  dx_nnz);
		}
	}
}

static struct var_type RECV_VAR_IRECV_REP = {
	VAR_RECV_VAR,
	irecv_init,
	NULL,			// deinit
	{
	 irecv_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL			// clear
	 }
};

const struct var_type *RECV_VAR_IRECV = &RECV_VAR_IRECV_REP;
