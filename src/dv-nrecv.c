#include "port.h"
#include <assert.h>
#include <string.h>
#include "frame.h"
#include "vars.h"

static void nrecv_init(struct design_var *v, const struct design *d,
		       void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	size_t n = frame_interval_count(design_frame(d));
	v->dim = n + 1;
	v->names = var_names_alloc("NRecv", strlen("NRecv"), n + 1);
}

static void nrecv_deinit(struct design_var *v)
{
	var_names_free(v->names);
}

static void nrecv_message_add(void *udata, struct frame *f,
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

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		size_t ix = frame_dyad_ix(f, isend, jrecv);
		design_update(d, ix, dx_data, &pat);
	}
}

static void nrecv_message_advance(void *udata, struct frame *f,
				  const struct message *msg, size_t intvl)
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

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	struct vpattern pat = vpattern_make(dx_index, 2);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		size_t ix1 = dyn_index + intvl;
		size_t ix0 = ix1 - 1;

		dx_index[0] = ix0;
		dx_index[1] = ix1;

		size_t ix = frame_dyad_ix(f, isend, jrecv);
		design_update(d, ix, dx_data, &pat);
	}
}

static struct var_type DYAD_VAR_NRECV_REP = {
	VAR_DYAD_VAR,
	nrecv_init,
	nrecv_deinit,
	{
	 nrecv_message_add,
	 nrecv_message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *DYAD_VAR_NRECV = &DYAD_VAR_NRECV_REP;
