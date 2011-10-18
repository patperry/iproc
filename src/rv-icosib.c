#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  v---- i
 *  k
 *  ^---- j
 */

static char *icosib_names[] = { "ICosib" };

static void icosib_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = icosib_names;
}

static void icosib_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim <=
			design_dvars_dim(frame_recv_design(f)));

	size_t dyn_index = v->dyn_index;
	size_t *imsg, i, n;

	double dx_data[1] = { +1.0 };
	ssize_t dx_index[1] = { dyn_index };
	size_t dx_nnz = 1;
	size_t dx_n = design_dvars_dim(frame_recv_design(f));
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	size_t isend = msg->from;
	size_t cojrecv = msg->from;

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;

		size_t krecv = msg->to[ito];

		frame_get_recv_messages(f, krecv, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;

			assert(msg1->to[ito] == msg1->to[ito]);
			if (msg1->from == msg->from
			    || msg1->from == msg->to[ito])
				continue;

			size_t jrecv = msg1->from;
			size_t coisend = msg1->from;

			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);

			const struct vector *dx =
			    frame_recv_dx(f, isend, jrecv);
			if (vector_item(dx, dyn_index) == 0.0) {
				frame_recv_update(f, isend, jrecv, &delta);
			}

			assert(coisend != cojrecv);
			assert(coisend != krecv);
			assert(cojrecv != krecv);

			const struct vector *codx =
			    frame_recv_dx(f, coisend, cojrecv);
			if (vector_item(codx, dyn_index) == 0.0) {
				frame_recv_update(f, coisend, cojrecv, &delta);
			}
		}
	}
}

static struct var_type RECV_VAR_ICOSIB_REP = {
	VAR_RECV_VAR,
	icosib_init,
	NULL,			// deinit
	{
	 icosib_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// send_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_ICOSIB = &RECV_VAR_ICOSIB_REP;
