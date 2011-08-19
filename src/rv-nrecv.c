#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void nrecv_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
}

static void nrecv_message_add(void *udata, struct frame *f,
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

	double dx_data = +1.0;
	ssize_t dx_index = dyn_index;
	ssize_t dx_nnz = 1;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(&dx_index, &dx_data, dx_nnz, dx_n);

	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t isend = msg->to[ito];

		frame_recv_update(f, isend, jrecv, &delta);
	}
}

static void nrecv_message_advance(void *udata, struct frame *f,
				  const struct message *msg, ssize_t intvl)
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

	double dx_data[2] = { -1.0, +1.0 };
	ssize_t dx_index[2] = { 0, 1 };
	ssize_t dx_nnz = 2;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t isend = msg->to[ito];
		ssize_t ix1 = dyn_index + intvl;
		ssize_t ix0 = ix1 - 1;

		dx_index[0] = ix0;
		dx_index[1] = ix1;
		
		frame_recv_update(f, isend, jrecv, &delta);
	}
}

static struct var_type RECV_VAR_NRECV_REP = {
	VAR_RECV_VAR,
	nrecv_init,
	NULL, // deinit
	NULL, // frame_init
	NULL, // frame_deinit
	{
		nrecv_message_add,
		nrecv_message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL,			// clear
	}
};

const struct var_type *RECV_VAR_NRECV = &RECV_VAR_NRECV_REP;
