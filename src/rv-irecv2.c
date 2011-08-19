#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void irecv2_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

	dv->dim = 1;
}

static void irecv2_message_add(void *udata, struct frame *f,
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
	ssize_t *imsg, i, n;	
	
	double dx_data[1] = { +1.0 };
	ssize_t dx_index[1] = { dyn_index };
	ssize_t dx_nnz = 1;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t isend = msg->to[ito];

		frame_get_send_messages(f, isend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			assert(msg1->from == isend);
			
			ssize_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				ssize_t isend1 = msg1->to[ito1];				
				const struct vector *dx = frame_recv_dx(f, isend1, jrecv);
				
				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, isend1, jrecv, &delta);
				}
			}
		}
	}
		
	frame_get_recv_messages(f, jrecv, &imsg, &n);
	for (i = 0; i < n; i++) {
		const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
		const struct message *msg1 = fmsg->message;
		ssize_t jrecv1 = msg1->from;
		
		for (ito = 0; ito < nto; ito++) {
			ssize_t isend = msg->to[ito];
			const struct vector *dx = frame_recv_dx(f, isend, jrecv1);
			
			if (vector_item(dx, dyn_index) == 0.0) {
				frame_recv_update(f, isend, jrecv1, &delta);
			}
		}
	}
}

static struct var_type RECV_VAR_IRECV2_REP = {
	VAR_RECV_VAR,
	irecv2_init,
	NULL, // deinit
	NULL, // frame_init
	NULL, // frame_deinit
	{
		irecv2_message_add,
		NULL,			// message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL,			// clear
	}
};

const struct var_type *RECV_VAR_IRECV2 = &RECV_VAR_IRECV2_REP;
