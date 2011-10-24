#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  i ----> k ----> j
 */

static char *isend2_names[] = { "ISend2" };

static void isend2_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = isend2_names;
}

static void isend2_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t dyn_index = v->dyn_index;
	size_t *imsg, i, n;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);

	size_t isend = msg->from;
	size_t cokrecv = msg->from;

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t ksend = msg->to[ito];

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->from)
					continue;

				size_t jrecv = msg1->to[ito1];

				assert(isend != jrecv);
				assert(isend != ksend);
				assert(jrecv != ksend);

				const struct vector *dx =
				    frame_recv_dx(f, isend, jrecv);

				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, isend, jrecv,
							  dx_data, &pat);
				}
			}
		}
	}

	frame_get_recv_messages(f, cokrecv, &imsg, &n);
	for (i = 0; i < n; i++) {
		const struct frame_message *fmsg =
		    frame_messages_item(f, imsg[i]);
		const struct message *msg1 = fmsg->message;

		if (msg1->from == msg->from)
			continue;

		size_t intvl1 = fmsg->interval;
		size_t ix = dyn_index + intvl1;
		size_t coisend = msg1->from;

		dx_index[0] = ix;

		for (ito = 0; ito < nto; ito++) {
			if (msg->from == msg->to[ito]
			    || msg1->from == msg->to[ito])
				continue;

			size_t cojrecv = msg->to[ito];

			assert(coisend != cojrecv);
			assert(coisend != cokrecv);
			assert(cojrecv != cokrecv);

			const struct vector *codx =
			    frame_recv_dx(f, coisend, cojrecv);

			if (vector_item(codx, dyn_index) == 0.0) {
				frame_recv_update(f, coisend, cojrecv,
						  dx_data, &pat);
			}
		}
	}
}

static struct var_type RECV_VAR_ISEND2_REP = {
	VAR_RECV_VAR,
	isend2_init,
	NULL,			// deinit
	{
	 isend2_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_ISEND2 = &RECV_VAR_ISEND2_REP;
