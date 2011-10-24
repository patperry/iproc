#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  i <---- k <---- j
 */

static char *irecv2_names[] = { "IRecv2" };

static void irecv2_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = irecv2_names;
}

static void irecv2_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim <= design_dvars_dim(frame_recv_design(f)));

	size_t dyn_index = v->dyn_index;
	double dx_data[1] = { 1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);
	size_t *imsg, i, n;

	size_t krecv = msg->from;
	size_t cojrecv = msg->from;

	size_t ito, nto = msg->nto;

	frame_get_recv_messages(f, krecv, &imsg, &n);
	for (i = 0; i < n; i++) {
		const struct frame_message *fmsg =
		    frame_messages_item(f, imsg[i]);
		const struct message *msg1 = fmsg->message;

		if (msg1->from == msg->from)
			continue;

		size_t jrecv = msg1->from;

		for (ito = 0; ito < nto; ito++) {
			if (msg->to[ito] == msg->from
			    || msg->to[ito] == msg1->from)
				continue;

			size_t isend = msg->to[ito];

			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);

			const struct vector *dx =
			    frame_recv_dx(f, isend, jrecv);

			if (vector_item(dx, dyn_index) == 0.0) {
				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
	}

	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t coksend = msg->to[ito];

		frame_get_send_messages(f, coksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->from)
					continue;

				size_t coisend = msg1->to[ito1];

				assert(coisend != cojrecv);
				assert(coisend != coksend);
				assert(cojrecv != coksend);

				const struct vector *dx =
				    frame_recv_dx(f, coisend, cojrecv);

				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, coisend, cojrecv,
							  dx_data, &pat);
				}
			}
		}
	}
}

static struct var_type RECV_VAR_IRECV2_REP = {
	VAR_RECV_VAR,
	irecv2_init,
	NULL,			// deinit
	{
	 irecv2_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_IRECV2 = &RECV_VAR_IRECV2_REP;
