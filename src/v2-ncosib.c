#include "port.h"
#include <assert.h>
#include <string.h>
#include "frame.h"
#include "vars.h"

/*
 *    I1
 *  v---- i
 *  h
 *  ^---- j
 *    I2
 */


static void ncosib_message_add(void *udata, struct history *h,
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

	const size_t *indx;
	size_t iz, nz;

	size_t isend = msg->from;
	size_t cojrecv = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t hrecv = msg->to[ito];
		history_get_recv_active(h, hrecv, &indx, &nz);
		const size_t *nmsg = history_recv_counts(h, hrecv);

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = indx[iz];
			size_t coisend = jrecv;

			for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t ix = iintvl * nintvl;
				size_t coix = iintvl;

				dx_data[0] = +(double)(*nmsg);

				if (hrecv != isend && hrecv != jrecv && isend != jrecv) {
					dx_index[0] = ix;
					design2_update(d, v, isend, jrecv, dx_data, dx_index, dx_nz);
				}

				if (!hrecv != coisend && hrecv != cojrecv && coisend != cojrecv) {
					dx_index[0] = coix;
					design2_update(d, v, coisend, cojrecv, dx_data, dx_index, dx_nz);
				}
			}
		}
	}
}


static void ncosib_message_advance(void *udata, struct history *h,
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

	size_t isend = msg->from;
	size_t cojrecv = msg->from;

	for (ito = 0; ito < nto; ito++) {
		size_t hrecv = msg->to[ito];
		history_get_recv_active(h, hrecv, &indx, &nz);
		const size_t *nmsg = history_recv_counts(h, hrecv);

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = indx[iz];
			size_t coisend = jrecv;

			for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
				if (*nmsg == 0)
					continue;

				size_t ix0 = (intvl - 1) + iintvl * nintvl;
				size_t ix1 = ix0 + 1;
				size_t coix0 = iintvl + (intvl - 1) * nintvl;
				size_t coix1 = coix0 + nintvl;

				dx_data[0] = -(double)(*nmsg);
				dx_data[1] = +(double)(*nmsg);

				if (hrecv != isend && hrecv != jrecv && isend != jrecv) {
					dx_index[0] = ix0;
					dx_index[1] = ix1;
					design2_update(d, v, isend, jrecv, dx_data, dx_index, dx_nz);
				}

				if (!hrecv != coisend && hrecv != cojrecv && coisend != cojrecv) {
					dx_index[0] = coix0;
					dx_index[1] = coix1;
					design2_update(d, v, coisend, cojrecv, dx_data, dx_index, dx_nz);
				}
			}
		}
	}
}


static struct history_callbacks ncosib_history_callbacks = {
	ncosib_message_add,
	ncosib_message_advance,
	NULL
};


static void ncosib_init(struct tvar2 *tv, struct design2 *d, va_list ap)
{
	(void)ap;		// unused;

	struct history *h = design2_history(d);
	size_t n = history_interval_count(h);

	tv->var.rank = 2;
	tv->var.dims[0] = n;
	tv->var.dims[1] = n;
	tv->udata = NULL;

	history_add_observer(h, tv, &ncosib_history_callbacks);
}


static void ncosib_deinit(struct tvar2 *tv, struct design2 *d)
{
	struct history *h = design2_history(d);
	history_remove_observer(h, tv);
}


static struct tvar2_type VAR2_NCOSIB_REP = {
	ncosib_init,
	ncosib_deinit,
};


const struct tvar2_type *VAR2_NCOSIB = &VAR2_NCOSIB_REP;

