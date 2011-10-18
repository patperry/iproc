#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *   /----> i
 *  k
 *   \----> j
 */

static char *isib_names[] = { "ISib" };

static void isib_init(struct design_var *v, const struct design *d,
		      void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = isib_names;
}

static void isib_message_add(void *udata, struct frame *f,
			     const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t ksend = msg->from;
	size_t dyn_index = v->dyn_index;
	size_t *imsg, i, n;

	double dx_data[1] = { +1.0 };
	ssize_t dx_index[1] = { dyn_index };
	size_t dx_nnz = 1;
	size_t dx_n = design_dvars_dim(frame_recv_design(f));
	struct svector delta = svector_make(dx_index, dx_data, dx_nnz, dx_n);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		size_t cojrecv = msg->to[ito];

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;
			assert(msg1->from == ksend);

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->to[ito])
					continue;

				size_t jrecv = msg1->to[ito1];

				assert(isend != jrecv);
				assert(isend != ksend);
				assert(jrecv != ksend);

				const struct vector *dx =
				    frame_recv_dx(f, isend, jrecv);
				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, isend, jrecv,
							  &delta);
				}
			}

			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->to[ito])
					continue;

				size_t coisend = msg1->to[ito1];

				assert(coisend != cojrecv);
				assert(coisend != ksend);
				assert(cojrecv != ksend);

				const struct vector *dx =
				    frame_recv_dx(f, coisend, cojrecv);
				if (vector_item(dx, dyn_index) == 0.0) {
					frame_recv_update(f, coisend, cojrecv,
							  &delta);
				}
			}
		}
	}
}

static struct var_type RECV_VAR_ISIB_REP = {
	VAR_RECV_VAR,
	isib_init,
	NULL,			// deinit
	{
	 isib_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// send_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_ISIB = &RECV_VAR_ISIB_REP;
