#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void nsend_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(dv);
	assert(d);
	assert(!params);

	size_t n = design_interval_count(d);
	dv->dim = n + 1;
	dv->names = var_names_alloc("NSend", strlen("NSend"), n + 1);
}

static void nsend_deinit(struct design_var *dv)
{
	var_names_free(dv->names);
}

static void nsend_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	struct frame_var *fv = udata;

	assert(fv);
	assert(f);
	assert(msg);
	assert(fv->design);
	assert(fv->design->dyn_index + fv->design->dim
	       <= design_dvars_dim(f->design));

	size_t isend = msg->from;
	size_t dyn_index = fv->design->dyn_index;

	double dx_data[1] = { +1.0 };
	ssize_t dx_index[1] = { dyn_index };
	size_t dx_nnz = 1;
	size_t dx_n = design_dvars_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;

		size_t jrecv = msg->to[ito];

		frame_recv_update(f, isend, jrecv, &delta);
	}
}

static void nsend_message_advance(void *udata, struct frame *f,
				  const struct message *msg, size_t intvl)
{
	struct frame_var *fv = udata;

	assert(fv);
	assert(f);
	assert(msg);
	assert(fv->design);
	assert(fv->design->dyn_index + fv->design->dim
	       <= design_dvars_dim(f->design));

	size_t isend = msg->from;
	size_t dyn_index = fv->design->dyn_index;

	double dx_data[2] = { -1.0, +1.0 };
	ssize_t dx_index[2] = { 0, 1 }; // values are unused
	size_t dx_nnz = 2;
	size_t dx_n = design_dvars_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;

		size_t jrecv = msg->to[ito];

		ssize_t ix1 = dyn_index + intvl;
		ssize_t ix0 = ix1 - 1;

		dx_index[0] = ix0;
		dx_index[1] = ix1;

		frame_recv_update(f, isend, jrecv, &delta);
	}
}

static struct var_type RECV_VAR_NSEND_REP = {
	VAR_RECV_VAR,
	nsend_init,
	nsend_deinit,
	NULL,			// frame_init
	NULL,			// frame_deinit
	{
	 nsend_message_add,
	 nsend_message_advance,
	 NULL,			// recv_update
	 NULL,			// send_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_NSEND = &RECV_VAR_NSEND_REP;
