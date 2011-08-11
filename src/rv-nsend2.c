#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

static void nsend2_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

	ssize_t n = vector_dim(design_intervals(d));
	ssize_t n1 = n + 1;
	dv->dim = n1 * n1;
}

static void handle_message_add(void *udata, struct frame *f,
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

	ssize_t isend = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;

	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t jrecv = msg->to[ito];

		ssize_t i, n = fv->design->dim;
		for (i = 0; i < n; i++) {
			frame_recv_update(f, isend, jrecv, dyn_index + i, 1.0);
		}
	}
}

static void handle_message_advance(void *udata, struct frame *f,
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
	ssize_t isend = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;

	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t jrecv = msg->to[ito];
		ssize_t *imsg, i, n;
		
		frame_get_send_messages(f, jrecv, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			ssize_t intvl1 = fmsg->interval;
			ssize_t ix = dyn_index + (intvl - 1) + intvl1 * (nintvl + 1);
			assert(msg1->from == jrecv);
			
			ssize_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				ssize_t jrecv1 = msg1->to[ito1];
				
				frame_recv_update(f, isend, jrecv1, ix, -1.0);
			}
		}

		frame_get_recv_messages(f, isend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg = frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			ssize_t intvl1 = fmsg->interval;
			ssize_t ix = dyn_index + intvl1 + (intvl - 1) * (nintvl + 1);
			
			ssize_t isend1 = msg1->from;
			frame_recv_update(f, isend1, jrecv, ix, -1.0);
		}
	}
}

static struct var_type RECV_VAR_NSEND2_REP = {
	VAR_RECV_VAR,
	nsend2_init,
	NULL,
	NULL,
	NULL,
	{
	 handle_message_add,
	 handle_message_advance,
	 NULL,			// recv_update
	 NULL,			// send_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_NSEND2 = &RECV_VAR_NSEND2_REP;
