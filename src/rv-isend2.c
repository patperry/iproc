#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  i ----> k ----> j
 */

static char *isend2_names[] = { "ISend2" };

static void isend2_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	(void)d; // unused
	(void)params; // unused;
	assert(dv);
	assert(d);
	assert(!params);

	dv->dim = 1;
	dv->names = isend2_names;
}


static void isend2_message_add(void *udata, struct frame *f,
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

	ssize_t dyn_index = fv->design->dyn_index;
	ssize_t *imsg, i, n;	
	
	double dx_data[1] = { +1.0 };
	ssize_t dx_index[1] = { dyn_index };
	ssize_t dx_nnz = 1;
	ssize_t dx_n = design_recv_dyn_dim(f->design);
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	ssize_t isend = msg->from;
	ssize_t cokrecv = isend;
	
	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t ksend = msg->to[ito];

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			
			if (msg1 == msg)
				continue;
			
			ssize_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				ssize_t jrecv = msg1->to[ito1];
				const struct vector *dx = frame_recv_dx(f, isend, jrecv);
				
				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, isend, jrecv, &delta);
				}
			}
		}
	}
		
	frame_get_recv_messages(f, cokrecv, &imsg, &n);
	for (i = 0; i < n; i++) {
		const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
		const struct message *msg1 = fmsg->message;
		
		if (msg1 == msg)
			continue;
		
		ssize_t intvl1 = fmsg->interval;
		ssize_t ix = dyn_index + intvl1;
		ssize_t coisend = msg1->from;
		
		dx_index[0] = ix;
		
		for (ito = 0; ito < nto; ito++) {
			ssize_t cojrecv = msg->to[ito];
			
			const struct vector *codx = frame_recv_dx(f, coisend, cojrecv);
			
			if (vector_item(codx, dyn_index) == 0.0) {
				frame_recv_update(f, coisend, cojrecv, &delta);
			}
		}
	}
}


static struct var_type RECV_VAR_ISEND2_REP = {
	VAR_RECV_VAR,
	isend2_init,
	NULL, // deinit
	NULL, // frame_init
	NULL, // frame_deinit
	{
		isend2_message_add,
		NULL,			// message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL,			// clear
	}
};

const struct var_type *RECV_VAR_ISEND2 = &RECV_VAR_ISEND2_REP;
