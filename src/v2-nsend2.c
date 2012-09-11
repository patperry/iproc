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

static void nsend2_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	if (!frame_interval_count(f))
		return;

	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	size_t dx_nz = 1;

	const struct history *h = frame_history(f);
	size_t iintvl, nintvl = frame_interval_count(f);
	size_t ito, nto = msg->nto;
	const size_t *nmsg;

	size_t iz, nz;
	const size_t *indx;
	size_t isend = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t hsend = msg->to[ito];
		history_get_send_active(h, hsend, &indx, &nz);
		nmsg = history_send_counts(h, hsend);

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = indx[iz];
			for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t ix = iintvl * nintvl;

				dx_data[0] = (double)(*nmsg);
				dx_index[0] = ix;

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				design2_update(d, v, isend, jrecv, dx_data, dx_index, dx_nz);
			}
		}
	}

	size_t cohrecv = isend;
	history_get_recv_active(h, cohrecv, &indx, &nz);
	nmsg = history_recv_counts(h, cohrecv);

	for (iz = 0; iz < nz; iz++) {
		size_t coisend = indx[iz];

		for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
			if (*nmsg == 0)
				continue;
			size_t coix = iintvl;

			dx_data[0] = +(double)(*nmsg);
			dx_index[0] = coix;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (cohrecv == coisend || cohrecv == cojrecv || coisend == cojrecv)
					continue;

				design2_update(d, v, coisend, cojrecv, dx_data, dx_index, dx_nz);
			}
		}
	}
}


static void nsend2_message_advance(void *udata, struct frame *f,
				   const struct message *msg, size_t intvl)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2];
	size_t dx_nz = 2;

	const struct history *h = frame_history(f);
	size_t iintvl, nintvl = frame_interval_count(f);
	size_t ito, nto = msg->nto;

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
			for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t ix0 = (intvl - 1) + iintvl * nintvl;
				size_t ix1 = ix0 + 1;

				dx_data[0] = -(double)(*nmsg);
				dx_data[1] = +(double)(*nmsg);
				dx_index[0] = ix0;
				dx_index[1] = ix1;

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				design2_update(d, v, isend, jrecv, dx_data, dx_index, dx_nz);
			}
		}
	}

	size_t cohrecv = isend;
	history_get_recv_active(h, cohrecv, &indx, &nz);
	nmsg = history_recv_counts(h, cohrecv);

	for (iz = 0; iz < nz; iz++) {
		size_t coisend = indx[iz];

		for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
			if (*nmsg == 0)
				continue;
			size_t coix0 = iintvl + (intvl - 1) * nintvl;
			size_t coix1 = coix0 + nintvl;

			dx_data[0] = -(double)(*nmsg);
			dx_data[1] = +(double)(*nmsg);
			dx_index[0] = coix0;
			dx_index[1] = coix1;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (cohrecv == coisend || cohrecv == cojrecv || coisend == cojrecv)
					continue;

				design2_update(d, v, coisend, cojrecv, dx_data, dx_index, dx_nz);
			}
		}
	}
}



static struct frame_callbacks nsend2_frame_callbacks = {
	nsend2_message_add,
	nsend2_message_advance,
	NULL
};


static void nsend2_init(struct tvar2 *tv, struct design2 *d, va_list ap)
{
	(void)ap;		// unused;

	struct frame *f = design2_frame(d);
	size_t n = frame_interval_count(f);

	tv->var.dim = n * n;
	tv->udata = NULL;

	frame_add_observer(f, tv, &nsend2_frame_callbacks);
}


static void nsend2_deinit(struct tvar2 *tv, struct design2 *d)
{
	struct frame *f = design2_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar2_type VAR2_NSEND2_REP = {
	nsend2_init,
	nsend2_deinit,
};


const struct tvar2_type *VAR2_NSEND2 = &VAR2_NSEND2_REP;



