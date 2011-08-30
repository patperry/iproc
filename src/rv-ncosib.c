#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *    I1
 *  v---- i
 *  k
 *  ^---- j
 *    I2
 */

static void ncosib_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	(void)d; // unused
	(void)params; // unused;
	assert(dv);
	assert(d);
	assert(!params);

	ssize_t n = vector_dim(design_intervals(d));
	ssize_t n1 = n + 1;
	dv->dim = n1 * n1;
	dv->names = var_names_alloc2("NCosib", strlen("NCosib"), n + 1, n + 1);
}

static void ncosib_deinit(struct design_var *dv)
{
	var_names_free(dv->names);
}

static void ncosib_message_add(void *udata, struct frame *f,
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

	const struct design *d = frame_design(f);
	const struct vector *intvls = design_intervals(d);
	ssize_t nintvl = vector_dim(intvls);

	ssize_t dyn_index = fv->design->dyn_index;
	ssize_t *imsg, i, n;	
	
	double dx_data[1] = { +1.0 };
	ssize_t dx_index[1] = { 0 };
	ssize_t dx_nnz = 1;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	ssize_t isend = msg->from;
	ssize_t cojrecv = msg->from;
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;
		
		ssize_t krecv = msg->to[ito];

		frame_get_recv_messages(f, krecv, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			
			assert(msg1->to[ito] == msg1->to[ito]);
			if (msg1->from == msg->from || msg1->from == msg->to[ito])
				continue;
			
			ssize_t intvl = fmsg->interval;
			ssize_t coix = dyn_index + intvl;			
			ssize_t ix = dyn_index + intvl * (nintvl + 1);			
			ssize_t jrecv = msg1->from;
			ssize_t coisend = msg1->from;			
			
			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);
			
			dx_index[0] = ix;
			frame_recv_update(f, isend, jrecv, &delta);
			
			assert(coisend != cojrecv);
			assert(coisend != krecv);
			assert(cojrecv != krecv);
			
			dx_index[0] = coix;
			frame_recv_update(f, coisend, cojrecv, &delta);
		}
	}
}

static void ncosib_message_advance(void *udata, struct frame *f,
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

	const struct design *d = frame_design(f);
	const struct vector *intvls = design_intervals(d);
	ssize_t nintvl = vector_dim(intvls);

	ssize_t dyn_index = fv->design->dyn_index;
	ssize_t ito, nto = msg->nto;
	ssize_t *imsg, i, n;

	double dx_data[2] = { -1.0, +1.0 };
	ssize_t dx_index[2] = { 0, 1 };
	ssize_t dx_nnz = 2;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);
	
	ssize_t isend = msg->from;
	ssize_t cojrecv = msg->from;
	
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;

		ssize_t krecv = msg->to[ito];
		
		frame_get_recv_messages(f, krecv, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			
			assert(msg1->to[ito] == msg1->to[ito]);
			if (msg1->from == msg->from || msg1->from == msg->to[ito])
				continue;
			
			ssize_t intvl1 = fmsg->interval;
			ssize_t ix0 = dyn_index + (intvl - 1) + intvl1 * (nintvl + 1);
			ssize_t ix1 = ix0 + 1;
			ssize_t coix0 = dyn_index + intvl1 + (intvl - 1) * (nintvl + 1);
			ssize_t coix1 = coix0 + (nintvl + 1);
			
			ssize_t jrecv = msg1->from;
			ssize_t coisend = jrecv;
			
			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);
			
			dx_index[0] = ix0;
			dx_index[1] = ix1;
			frame_recv_update(f, isend, jrecv, &delta);
			
			assert(coisend != cojrecv);
			assert(coisend != krecv);
			assert(cojrecv != krecv);
			
			dx_index[0] = coix0;
			dx_index[1] = coix1;
			frame_recv_update(f, coisend, cojrecv, &delta);
		}
	}
}

static struct var_type RECV_VAR_NCOSIB_REP = {
	VAR_RECV_VAR,
	ncosib_init,
	ncosib_deinit,
	NULL, // frame_init
	NULL, // frame_deinit
	{
		ncosib_message_add,
		ncosib_message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL,			// clear
	}
};

const struct var_type *RECV_VAR_NCOSIB = &RECV_VAR_NCOSIB_REP;
