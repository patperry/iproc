#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static char *irecv_names[] = { "IRecv" };

static void irecv_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	(void)d; // unused
	(void)params; // unused;
	assert(dv);
	assert(d);
	assert(!params);

	dv->dim = 1;
	dv->names = irecv_names;
}

static void irecv_message_add(void *udata, struct frame *f,
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

	ssize_t jrecv = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;
	ssize_t dx_index = dyn_index;
	double dx_data = 1.0;
	ssize_t dx_nnz = 1;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(&dx_index, &dx_data, dx_nnz, dx_n);

	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t isend = msg->to[ito];
		const struct vector *dx = frame_recv_dx(f, isend, jrecv);
		if (vector_item(dx, dyn_index) == 0.0) {
			frame_recv_update(f, isend, jrecv, &delta);
		}
	}
}

static struct var_type RECV_VAR_IRECV_REP = {
	VAR_RECV_VAR,
	irecv_init,
	NULL, // deinit
	NULL, // frame_init
	NULL, // frame_deinit
	{
		irecv_message_add,
		NULL,			// message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL			// clear
	}
};

const struct var_type *RECV_VAR_IRECV = &RECV_VAR_IRECV_REP;
