#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *   /----> i
 *  k
 *   \----> j
 */

static void isib_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	(void)d; // unused
	(void)params; // unused;
	assert(dv);
	assert(d);
	assert(!params);

	dv->dim = 1;
}

static void isib_message_add(void *udata, struct frame *f,
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
	ssize_t ksend = msg->from;
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
		ssize_t cojrecv = isend;

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			ssize_t intvl = fmsg->interval;
			ssize_t ix = dyn_index + intvl * (nintvl + 1);
			ssize_t coix = dyn_index + intvl;
			assert(msg1->from == ksend);
			
			dx_index[0] = ix;
			
			ssize_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				ssize_t jrecv = msg1->to[ito1];
				const struct vector *dx = frame_recv_dx(f, isend, jrecv);
				
				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, isend, jrecv, &delta);
				}
			}
			
			dx_index[0] = coix;
			
			for (ito1 = 0; ito1 < nto1; ito1++) {
				ssize_t coisend = msg1->to[ito1];
				const struct vector *dx = frame_recv_dx(f, coisend, cojrecv);
				
				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, coisend, cojrecv, &delta);
				}
			}
		}
	}
}


static struct var_type RECV_VAR_ISIB_REP = {
	VAR_RECV_VAR,
	isib_init,
	NULL, // deinit
	NULL, // frame_init
	NULL, // frame_deinit
	{
		isib_message_add,
		NULL,                   // message_advance,
		NULL,			// recv_update
		NULL,			// send_update
		NULL,			// clear
	}
};

const struct var_type *RECV_VAR_ISIB = &RECV_VAR_ISIB_REP;
