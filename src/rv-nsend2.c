#include "port.h"
#include <assert.h>
#include <string.h>
#include "frame.h"
#include "vars.h"

/*
 *     I1      I2
 *  i ----> h ----> j
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

	const struct history *h = frame_history(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	size_t dyn_index = v->dyn_index;
	size_t ito, nto = msg->nto;
	const size_t *nmsg;
	size_t intvl1;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	struct vpattern pat = vpattern_make(dx_index, 1);

	size_t iz, nz;
	const size_t *indx;
	size_t isend = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t hsend = msg->to[ito];
		history_get_send_active(h, hsend, &indx, &nz);
		nmsg = history_send_counts(h, hsend);

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = indx[iz];
			for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t ix = dyn_index + intvl1 * nintvl1;

				dx_data[0] = (double)(*nmsg);
				dx_index[0] = ix;

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
	}

	size_t cohrecv = isend;
	history_get_recv_active(h, cohrecv, &indx, &nz);
	nmsg = history_recv_counts(h, cohrecv);

	for (iz = 0; iz < nz; iz++) {
		size_t coisend = indx[iz];

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;
			size_t coix = dyn_index + intvl1;

			dx_data[0] = +(double)(*nmsg);
			dx_index[0] = coix;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (cohrecv == coisend || cohrecv == cojrecv || coisend == cojrecv)
					continue;

				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
			}
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

	const struct history *h = frame_history(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	size_t dyn_index = v->dyn_index;
	size_t ito, nto = msg->nto;
	size_t intvl1;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	struct vpattern pat = vpattern_make(dx_index, 2);

	size_t iz, nz;
	const size_t *indx;
	const size_t *nmsg;
	size_t isend = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t hsend = msg->to[ito];
		history_get_send_active(h, hsend, &indx, &nz);
		nmsg = history_send_counts(h, hsend);

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = indx[iz];
			for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t ix0 = dyn_index + (intvl - 1) + intvl1 * nintvl1;
				size_t ix1 = ix0 + 1;

				dx_data[0] = -(double)(*nmsg);
				dx_data[1] = +(double)(*nmsg);
				dx_index[0] = ix0;
				dx_index[1] = ix1;

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
	}

	size_t cohrecv = isend;
	history_get_recv_active(h, cohrecv, &indx, &nz);
	nmsg = history_recv_counts(h, cohrecv);

	for (iz = 0; iz < nz; iz++) {
		size_t coisend = indx[iz];

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;
			size_t coix0 = dyn_index + intvl1 + (intvl - 1) * nintvl1;
			size_t coix1 = coix0 + nintvl1;

			dx_data[0] = -(double)(*nmsg);
			dx_data[1] = +(double)(*nmsg);
			dx_index[0] = coix0;
			dx_index[1] = coix1;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (cohrecv == coisend || cohrecv == cojrecv || coisend == cojrecv)
					continue;

				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
			}
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
