#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static char *isend_names[] = { "ISend" };

static void isend_init(struct design_var *v, const struct design *d,
		       void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = isend_names;
}

static void isend_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t isend = msg->from;
	size_t dyn_index = v->dyn_index;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;

		size_t jrecv = msg->to[ito];
		const double *dx = frame_recv_dx(f, isend, jrecv);
		if (dx[dyn_index] == 0.0) {
			frame_recv_update(f, isend, jrecv, dx_data, &pat);
		}
	}
}

static struct var_type RECV_VAR_ISEND_REP = {
	VAR_RECV_VAR,
	isend_init,
	NULL,			// deinit
	{
	 isend_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL			// clear
	 }
};

const struct var_type *RECV_VAR_ISEND = &RECV_VAR_ISEND_REP;
