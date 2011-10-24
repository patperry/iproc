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

static void ncosib_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	const struct frame *f = design_frame(d);
	size_t n = frame_interval_count(f);
	size_t n1 = n + 1;
	v->dim = n1 * n1;
	v->names = var_names_alloc2("NCosib", strlen("NCosib"), n + 1, n + 1);
}

static void ncosib_deinit(struct design_var *v)
{
	var_names_free(v->names);
}

static void ncosib_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t nintvl = frame_interval_count(f);

	size_t dyn_index = v->dyn_index;
	size_t *imsg, i, n;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	struct vpattern pat = vpattern_make(dx_index, 1);

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

			size_t intvl = fmsg->interval;
			size_t coix = dyn_index + intvl;
			size_t ix = dyn_index + intvl * (nintvl + 1);
			size_t jrecv = msg1->from;
			size_t coisend = msg1->from;

			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);

			dx_index[0] = ix;
			frame_recv_update(f, isend, jrecv, dx_data, &pat);

			assert(coisend != cojrecv);
			assert(coisend != krecv);
			assert(cojrecv != krecv);

			dx_index[0] = coix;
			frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
		}
	}
}

static void ncosib_message_advance(void *udata, struct frame *f,
				   const struct message *msg, size_t intvl)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t nintvl = frame_interval_count(f);
	size_t dyn_index = v->dyn_index;
	size_t ito, nto = msg->nto;
	size_t *imsg, i, n;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	struct vpattern pat = vpattern_make(dx_index, 2);

	size_t isend = msg->from;
	size_t cojrecv = msg->from;

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

			size_t intvl1 = fmsg->interval;
			size_t ix0 =
			    dyn_index + (intvl - 1) + intvl1 * (nintvl + 1);
			size_t ix1 = ix0 + 1;
			size_t coix0 =
			    dyn_index + intvl1 + (intvl - 1) * (nintvl + 1);
			size_t coix1 = coix0 + (nintvl + 1);

			size_t jrecv = msg1->from;
			size_t coisend = jrecv;

			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);

			dx_index[0] = ix0;
			dx_index[1] = ix1;
			frame_recv_update(f, isend, jrecv, dx_data, &pat);

			assert(coisend != cojrecv);
			assert(coisend != krecv);
			assert(cojrecv != krecv);

			dx_index[0] = coix0;
			dx_index[1] = coix1;
			frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
		}
	}
}

static struct var_type RECV_VAR_NCOSIB_REP = {
	VAR_RECV_VAR,
	ncosib_init,
	ncosib_deinit,
	{
	 ncosib_message_add,
	 ncosib_message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_NCOSIB = &RECV_VAR_NCOSIB_REP;
