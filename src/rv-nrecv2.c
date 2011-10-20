#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *      I1      I2
 *  i <---- k <---- j
 *
 */

static void nrecv2_init(struct design_var *v, const struct design *d,
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
	v->names = var_names_alloc2("NRecv2", strlen("NRecv2"), n + 1, n + 1);
}

static void nrecv2_deinit(struct design_var *v)
{
	var_names_free(v->names);
}

static void nrecv2_message_add(void *udata, struct frame *f,
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
	double dx_data[1] = { 1.0 };
	size_t dx_index[1] = { dyn_index };
	size_t dx_nnz = 1;
	size_t *imsg, i, n;

	size_t krecv = msg->from;
	size_t cojrecv = krecv;

	size_t ito, nto = msg->nto;

	frame_get_recv_messages(f, krecv, &imsg, &n);
	for (i = 0; i < n; i++) {
		const struct frame_message *fmsg =
		    frame_messages_item(f, imsg[i]);
		const struct message *msg1 = fmsg->message;

		if (msg1->from == msg->from)
			continue;

		size_t intvl2 = fmsg->interval;
		size_t ix = dyn_index + intvl2 * (nintvl + 1);
		size_t jrecv = msg1->from;

		dx_index[0] = ix;

		for (ito = 0; ito < nto; ito++) {
			if (msg->to[ito] == msg->from
			    || msg->to[ito] == msg1->from)
				continue;

			size_t isend = msg->to[ito];

			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);

			frame_recv_update(f, isend, jrecv, dx_data, dx_index,
					  dx_nnz);
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

			size_t intvl1 = fmsg->interval;
			size_t coix = dyn_index + intvl1;

			dx_index[0] = coix;	// affects delta

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->from)
					continue;

				size_t coisend = msg1->to[ito1];

				assert(coisend != cojrecv);
				assert(coisend != coksend);
				assert(cojrecv != coksend);

				frame_recv_update(f, coisend, cojrecv, dx_data,
						  dx_index, dx_nnz);
			}
		}
	}

}

static void nrecv2_message_advance(void *udata, struct frame *f,
				   const struct message *msg, size_t intvl)
{
	struct design_var *v = udata;

	assert(f);
	assert(msg);
	assert(v);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t nintvl = frame_interval_count(f);
	size_t dyn_index = v->dyn_index;

	size_t ito, nto = msg->nto;
	size_t *imsg, i, n;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	size_t dx_nnz = 2;

	size_t krecv = msg->from;
	size_t cojrecv = krecv;

	frame_get_recv_messages(f, krecv, &imsg, &n);
	for (i = 0; i < n; i++) {
		const struct frame_message *fmsg =
		    frame_messages_item(f, imsg[i]);
		const struct message *msg1 = fmsg->message;

		if (msg1->from == msg->from)
			continue;

		size_t intvl2 = fmsg->interval;
		size_t ix0 = dyn_index + (intvl - 1) + intvl2 * (nintvl + 1);
		size_t ix1 = ix0 + 1;
		size_t jrecv = msg1->from;

		dx_index[0] = ix0;	// affects delta 
		dx_index[1] = ix1;	//

		for (ito = 0; ito < nto; ito++) {
			if (msg->to[ito] == msg->from
			    || msg->to[ito] == msg1->from)
				continue;

			size_t isend = msg->to[ito];

			assert(isend != jrecv);
			assert(isend != krecv);
			assert(jrecv != krecv);

			frame_recv_update(f, isend, jrecv, dx_data, dx_index,
					  dx_nnz);
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

			size_t intvl1 = fmsg->interval;
			size_t coix0 =
			    dyn_index + intvl1 + (intvl - 1) * (nintvl + 1);
			size_t coix1 = coix0 + (nintvl + 1);

			dx_index[0] = coix0;	// affects delta 
			dx_index[1] = coix1;	//

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->from)
					continue;

				size_t coisend = msg1->to[ito1];

				assert(coisend != cojrecv);
				assert(coisend != coksend);
				assert(cojrecv != coksend);

				frame_recv_update(f, coisend, cojrecv, dx_data,
						  dx_index, dx_nnz);
			}
		}
	}

}

static struct var_type RECV_VAR_NRECV2_REP = {
	VAR_RECV_VAR,
	nrecv2_init,
	nrecv2_deinit,
	{
	 nrecv2_message_add,
	 nrecv2_message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_NRECV2 = &RECV_VAR_NRECV2_REP;
