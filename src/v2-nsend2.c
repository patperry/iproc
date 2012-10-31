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

static void nsend2_message_add(void *udata, struct history *h,
			       const struct message *msg)
{
	if (!history_interval_count(h))
		return;

	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = v->design;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	size_t dx_nz = 1;

	size_t iintvl, nintvl = history_interval_count(h);
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


static void nsend2_message_advance(void *udata, struct history *h,
				   const struct message *msg, size_t intvl)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = v->design;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2];
	size_t dx_nz = 2;

	size_t iintvl, nintvl = history_interval_count(h);
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



static struct history_callbacks nsend2_history_callbacks = {
	nsend2_message_add,
	nsend2_message_advance,
	NULL
};


static void nsend2_init(struct tvar2 *tv, struct design2 *d, va_list ap)
{
	(void)ap;		// unused;

	struct frame *f = design2_frame(d);
	struct history *h = frame_history(f);
	size_t n = frame_interval_count(f);

	tv->var.rank = 2;
	tv->var.dims[0] = n;
	tv->var.dims[1] = n;
	tv->udata = NULL;

	history_add_observer(h, tv, &nsend2_history_callbacks);
}


static void nsend2_deinit(struct tvar2 *tv, struct design2 *d)
{
	struct frame *f = design2_frame(d);
	struct history *h = frame_history(f);
	history_remove_observer(h, tv);
}


static struct tvar2_type VAR2_NSEND2_REP = {
	nsend2_init,
	nsend2_deinit,
};


const struct tvar2_type *VAR2_NSEND2 = &VAR2_NSEND2_REP;



