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
	       <= design_dvars_dim(frame_dyad_design(f)));

	struct design *d = frame_dyad_design(f);
	size_t jrecv = msg->from;
	size_t dyn_index = v->dyn_index;
	double dx_data[1] = { 1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		size_t ix = frame_dyad_ix(f, isend, jrecv);
		const double *dx = design_dx(d, ix);
		if (!dx || dx[dyn_index] == 0.0) {
			design_update(d, ix, dx_data, &pat);
		}
	}
}

static struct var_type DYAD_VAR_IRECV_REP = {
	VAR_DYAD_VAR,
	irecv_init,
	NULL,			// deinit
	{
	 irecv_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL			// clear
	 }
};

const struct var_type *DYAD_VAR_IRECV = &DYAD_VAR_IRECV_REP;
