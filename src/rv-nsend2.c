#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *     I1      I2
 *  i ----> k ----> j
 *
 */

static void nsend2_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	size_t n = frame_interval_count(design_frame(d));
	size_t n1 = n + 1;
	v->dim = n1 * n1;
	v->names = var_names_alloc2("NSend2", strlen("NSend2"), n + 1, n + 1);
}

static void nsend2_deinit(struct design_var *v)
{
	var_names_free(v->names);
}

static void nsend2_message_add(void *udata, struct frame *f,
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
	size_t dx_nnz = 1;

	size_t isend = msg->from;
	size_t cokrecv = isend;

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

			size_t intvl1 = fmsg->interval;
			ssize_t ix = dyn_index + intvl1 * (nintvl + 1);

			dx_index[0] = ix;

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->from)
					continue;

				size_t jrecv = msg1->to[ito1];

				assert(isend != jrecv);
				assert(isend != ksend);
				assert(jrecv != ksend);

				frame_recv_update(f, isend, jrecv, dx_data,
						  dx_index, dx_nnz);
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
		ssize_t ix = dyn_index + intvl1;
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

			frame_recv_update(f, coisend, cojrecv, dx_data,
					  dx_index, dx_nnz);
		}
	}
}

static void nsend2_message_advance(void *udata, struct frame *f,
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
	size_t dx_nnz = 2;

	size_t isend = msg->from;
	size_t cokrecv = isend;

	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t ksend = msg->to[ito];

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;

			size_t intvl1 = fmsg->interval;
			ssize_t ix0 =
			    dyn_index + (intvl - 1) + intvl1 * (nintvl + 1);
			ssize_t ix1 = ix0 + 1;

			dx_index[0] = ix0;
			dx_index[1] = ix1;

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->from)
					continue;

				size_t jrecv = msg1->to[ito1];

				assert(isend != jrecv);
				assert(isend != ksend);
				assert(jrecv != ksend);

				frame_recv_update(f, isend, jrecv, dx_data,
						  dx_index, dx_nnz);
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
		ssize_t ix0 = dyn_index + intvl1 + (intvl - 1) * (nintvl + 1);
		ssize_t ix1 = ix0 + (nintvl + 1);
		size_t coisend = msg1->from;

		dx_index[0] = ix0;
		dx_index[1] = ix1;

		for (ito = 0; ito < nto; ito++) {
			if (msg->from == msg->to[ito]
			    || msg1->from == msg->to[ito])
				continue;

			size_t cojrecv = msg->to[ito];

			assert(coisend != cojrecv);
			assert(coisend != cokrecv);
			assert(cojrecv != cokrecv);

			frame_recv_update(f, coisend, cojrecv, dx_data,
					  dx_index, dx_nnz);
		}
	}
}

static struct var_type RECV_VAR_NSEND2_REP = {
	VAR_RECV_VAR,
	nsend2_init,
	nsend2_deinit,
	{
	 nsend2_message_add,
	 nsend2_message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_NSEND2 = &RECV_VAR_NSEND2_REP;
