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


static void nsib_message_add(void *udata, struct history *h,
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

		for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix = iintvl * nintvl;
			size_t coix = iintvl;

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

static void nsib_message_advance(void *udata, struct history *h,
				 const struct message *msg, size_t intvl)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = v->design;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2];
	size_t dx_nz = 2;

	size_t iintvl, nintvl = history_interval_count(h);
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

		for (iintvl = 0; iintvl < nintvl; iintvl++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix0 = (intvl - 1) + iintvl * nintvl;
			size_t ix1 = ix0 + 1;
			size_t coix0 = iintvl + (intvl - 1) * nintvl;
			size_t coix1 = coix0 + nintvl;

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


static struct history_callbacks nsib_history_callbacks = {
	nsib_message_add,
	nsib_message_advance,
	NULL
};


static void nsib_init(struct tvar2 *tv, const char *name, struct history *h, va_list ap)
{
	(void)ap;		// unused;

	size_t n = history_interval_count(h);
	size_t rank = 2;
	size_t dims[2] = { n, n };
	var_meta_init(&tv->var.meta, name, VAR_TYPE_TVAR, dims, rank);
	tv->udata = NULL;

	history_add_observer(h, tv, &nsib_history_callbacks);
}


static void nsib_deinit(struct tvar2 *tv, struct history *h)
{
	history_remove_observer(h, tv);
	var_meta_deinit(&tv->var.meta);
}


static struct tvar2_type VAR2_NSIB_REP = {
	nsib_init,
	nsib_deinit,
};


const struct tvar2_type *VAR2_NSIB = &VAR2_NSIB_REP;


