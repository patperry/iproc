#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *     I1
 *   /----> i
 *  k
 *   \----> j
 *     I2
 */

static void nsib_init(struct design_var *v, const struct design *d,
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
	v->names = var_names_alloc2("NSib", strlen("NSib"), n + 1, n + 1);
}

static void nsib_deinit(struct design_var *v)
{
	var_names_free(v->names);
}

static void nsib_message_add(void *udata, struct frame *f,
			     const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t nintvl = frame_interval_count(f);
	size_t ksend = msg->from;
	size_t dyn_index = v->dyn_index;
	size_t *imsg, i, n;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	struct vpattern pat = vpattern_make(dx_index, 1);

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		size_t cojrecv = isend;

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;

			size_t intvl = fmsg->interval;
			ssize_t ix = dyn_index + intvl * (nintvl + 1);
			ssize_t coix = dyn_index + intvl;
			assert(msg1->from == ksend);

			dx_index[0] = ix;

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->to[ito])
					continue;

				size_t jrecv = msg1->to[ito1];

				assert(isend != jrecv);
				assert(isend != ksend);
				assert(jrecv != ksend);

				frame_recv_update(f, isend, jrecv, dx_data,
						  &pat);
			}

			dx_index[0] = coix;

			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->to[ito])
					continue;

				size_t coisend = msg1->to[ito1];

				assert(coisend != cojrecv);
				assert(coisend != ksend);
				assert(cojrecv != ksend);

				frame_recv_update(f, coisend, cojrecv, dx_data,
						  &pat);
			}
		}
	}
}

static void nsib_message_advance(void *udata, struct frame *f,
				 const struct message *msg, size_t intvl)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t nintvl = frame_interval_count(f);
	size_t ksend = msg->from;
	size_t dyn_index = v->dyn_index;

	size_t ito, nto = msg->nto;
	size_t *imsg, i, n;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	struct vpattern pat = vpattern_make(dx_index, 2);

	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;

		size_t isend = msg->to[ito];
		size_t cojrecv = isend;

		frame_get_send_messages(f, ksend, &imsg, &n);
		for (i = 0; i < n; i++) {
			const struct frame_message *fmsg =
			    frame_messages_item(f, imsg[i]);
			const struct message *msg1 = fmsg->message;

			size_t intvl1 = fmsg->interval;
			ssize_t ix0 =
			    dyn_index + (intvl - 1) + intvl1 * (nintvl + 1);
			ssize_t ix1 = ix0 + 1;
			ssize_t coix0 =
			    dyn_index + intvl1 + (intvl - 1) * (nintvl + 1);
			ssize_t coix1 = coix0 + (nintvl + 1);
			assert(msg1->from == ksend);

			dx_index[0] = ix0;
			dx_index[1] = ix1;

			size_t ito1, nto1 = msg1->nto;
			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->to[ito])
					continue;

				size_t jrecv = msg1->to[ito1];

				assert(isend != jrecv);
				assert(isend != ksend);
				assert(jrecv != ksend);

				frame_recv_update(f, isend, jrecv, dx_data,
						  &pat);
			}

			dx_index[0] = coix0;
			dx_index[1] = coix1;

			for (ito1 = 0; ito1 < nto1; ito1++) {
				if (msg1->to[ito1] == msg1->from
				    || msg1->to[ito1] == msg->to[ito])
					continue;

				size_t coisend = msg1->to[ito1];

				assert(coisend != cojrecv);
				assert(coisend != ksend);
				assert(cojrecv != ksend);

				frame_recv_update(f, coisend, cojrecv, dx_data,
						  &pat);
			}
		}
	}

}

static struct var_type RECV_VAR_NSIB_REP = {
	VAR_RECV_VAR,
	nsib_init,
	nsib_deinit,
	{
	 nsib_message_add,
	 nsib_message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_NSIB = &RECV_VAR_NSIB_REP;
