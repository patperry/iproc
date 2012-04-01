#include "port.h"
#include <assert.h>
#include <string.h>
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

	const struct history *h = frame_history(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	size_t dyn_index = v->dyn_index;
	size_t ito, nto = msg->nto;

	double dx_data[1] = { 1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);


	const size_t *indx;
	size_t iz, nz;
	const size_t *nmsg;
	size_t intvl1;

	size_t hrecv = msg->from;
	history_get_recv_active(h, hrecv, &indx, &nz);
	nmsg = history_recv_counts(h, hrecv);

	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = indx[iz];

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix = dyn_index + intvl1 * nintvl1;

			dx_data[0] = +(double)(*nmsg);
			dx_index[0] = ix;

			for (ito = 0; ito < nto; ito++) {
				size_t isend = msg->to[ito];

				if (hrecv == isend || hrecv == jrecv || isend == jrecv)
					continue;

				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
	}

	size_t cojrecv = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t cohsend = msg->to[ito];
		history_get_send_active(h, cohsend, &indx, &nz);
		nmsg = history_send_counts(h, cohsend);

		for (iz = 0; iz < nz; iz++) {
			size_t coisend = indx[iz];

			for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t coix = dyn_index + intvl1;

				dx_data[0] = +(double)(*nmsg);
				dx_index[0] = coix;

				if (cohsend == coisend || cohsend == cojrecv || coisend == cojrecv)
					continue;

				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
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

	const struct history *h = frame_history(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	size_t dyn_index = v->dyn_index;

	size_t ito, nto = msg->nto;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	struct vpattern pat = vpattern_make(dx_index, 2);

	size_t iz, nz;
	const size_t *indx;
	const size_t *nmsg;
	size_t intvl1;

	size_t hrecv = msg->from;
	history_get_recv_active(h, hrecv, &indx, &nz);
	nmsg = history_recv_counts(h, hrecv);

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


			for (ito = 0; ito < nto; ito++) {
				size_t isend = msg->to[ito];

				if (hrecv == isend || hrecv == jrecv || isend == jrecv)
					continue;

				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
	}

	size_t cojrecv = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t cohsend = msg->to[ito];
		history_get_send_active(h, cohsend, &indx, &nz);
		nmsg = history_send_counts(h, cohsend);

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

				if (cohsend == coisend || cohsend == cojrecv || coisend == cojrecv)
					continue;

				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
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
