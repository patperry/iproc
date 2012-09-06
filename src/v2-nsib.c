#include "port.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "frame.h"
#include "vars.h"

/*
 *     I1
 *   /----> i
 *  h
 *   \----> j
 *     I2
 */


static void nsib_message_add(void *udata, struct frame *f,
			     const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	size_t dx_nz = 1;

	const struct history *h = frame_history(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	size_t hsend = msg->from;
	size_t ito, nto = msg->nto;

	size_t iz, nz;
	const size_t *indx;
	const size_t *nmsg;

	history_get_send_active(h, hsend, &indx, &nz);
	nmsg = history_send_counts(h, hsend);

	//printf("add { %zu -> [", msg->from);
	//for (ito = 0; ito < nto; ito++) {
	//	printf(" %zu", msg->to[ito]);
	//}
	//printf(" ] }\n");

	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = indx[iz];
		size_t coisend = jrecv;
		size_t intvl1;

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix = intvl1 * nintvl1;
			size_t coix = intvl1;

			dx_data[0] = +(double)(*nmsg);
			dx_index[0] = ix;

			for (ito = 0; ito < nto; ito++) {
				size_t isend = msg->to[ito];

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				design2_update(d, v, isend, jrecv, dx_data, dx_index, dx_nz);
				//printf("    dx[%zu, %zu](%zu, %zu) += %g\n",
				//       isend, jrecv, 0, intvl1, dx_data[0]);
			}

			dx_index[0] = coix;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (hsend == coisend || hsend == cojrecv || coisend == cojrecv)
					continue;

				design2_update(d, v, coisend, cojrecv, dx_data, dx_index, dx_nz);
				//printf("    dx[%zu, %zu](%zu, %zu) += %g\n",
				//       coisend, cojrecv, intvl1, 0, dx_data[0]);
			}

		}
	}


}

static void nsib_message_advance(void *udata, struct frame *f,
				 const struct message *msg, size_t intvl)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2];
	size_t dx_nz = 2;

	const struct history *h = frame_history(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	size_t hsend = msg->from;
	size_t ito, nto = msg->nto;

	size_t iz, nz;
	const size_t *indx;
	const size_t *nmsg;

	history_get_send_active(h, hsend, &indx, &nz);
	nmsg = history_send_counts(h, hsend);

	//printf("advance { %zu -> [", msg->from);
	//for (ito = 0; ito < nto; ito++) {
	//	printf(" %zu", msg->to[ito]);
	//}
	//printf(" ] } to %zu\n", intvl);


	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = indx[iz];
		size_t coisend = jrecv;
		size_t intvl1;

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix0 = (intvl - 1) + intvl1 * nintvl1;
			size_t ix1 = ix0 + 1;
			size_t coix0 = intvl1 + (intvl - 1) * nintvl1;
			size_t coix1 = coix0 + nintvl1;

			dx_data[0] = -(double)(*nmsg);
			dx_data[1] = +(double)(*nmsg);
			dx_index[0] = ix0;
			dx_index[1] = ix1;

			for (ito = 0; ito < nto; ito++) {
				size_t isend = msg->to[ito];

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				design2_update(d, v, isend, jrecv, dx_data, dx_index, dx_nz);
				//printf("    dx[%zu, %zu](%zu, %zu) += %g\n",
				//       isend, jrecv, intvl - 1, intvl1, dx_data[0]);
				//printf("    dx[%zu, %zu](%zu, %zu) += %g\n",
				//       isend, jrecv, intvl, intvl1, dx_data[1]);
			}

			dx_index[0] = coix0;
			dx_index[1] = coix1;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (hsend == coisend || hsend == cojrecv || coisend == cojrecv)
					continue;

				design2_update(d, v, coisend, cojrecv, dx_data, dx_index, dx_nz);
				//printf("    dx[%zu, %zu](%zu, %zu) += %g\n",
				//       coisend, cojrecv, intvl1, intvl - 1, dx_data[0]);
				//printf("    dx[%zu, %zu](%zu, %zu) += %g\n",
				//       coisend, cojrecv, intvl1, intvl, dx_data[1]);
			}

		}
	}

}


static struct frame_callbacks nsib_frame_callbacks = {
	nsib_message_add,
	nsib_message_advance,
	NULL
};


static void nsib_init(struct tvar2 *tv, const struct design2 *d, va_list ap)
{
	(void)ap;		// unused;

	struct frame *f = design2_frame(d);
	size_t n = frame_interval_count(f);

	tv->var.dim = (n + 1) * (n + 1);
	tv->udata = NULL;

	frame_add_observer(f, tv, &nsib_frame_callbacks);
}


static void nsib_deinit(struct tvar2 *tv, const struct design2 *d)
{
	struct frame *f = design2_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar2_type DYAD_VAR_NSIB_REP = {
	nsib_init,
	nsib_deinit,
};


const struct tvar2_type *DYAD_VAR_NSIB = &DYAD_VAR_NSIB_REP;


