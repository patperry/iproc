#include "port.h"
#include <assert.h>
#include <stdio.h>

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
	size_t nintvl1 = nintvl + 1;
	size_t hsend = msg->from;
	size_t dyn_index = v->dyn_index;

	size_t ito, nto = msg->nto;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0, };
	struct vpattern pat = vpattern_make(dx_index, 1);
	const struct frame_actor *fa = &f->senders[hsend];
	size_t iz, nz = fa->active.nz;
	const size_t *nmsg;

	//printf("add { %zd -> [", msg->from);
	//for (ito = 0; ito < nto; ito++) {
	//	printf(" %zd", msg->to[ito]);
	//}
	//printf(" ] }\n");

	for (iz = 0, nmsg = fa->nmsg; iz < nz; iz++) {
		size_t jrecv = fa->active.indx[iz];
		size_t coisend = jrecv;
		size_t intvl1;

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix = dyn_index + intvl1 * nintvl1;
			size_t coix = dyn_index + intvl1;

			dx_data[0] = +(double)(*nmsg);
			dx_index[0] = ix;

			for (ito = 0; ito < nto; ito++) {
				size_t isend = msg->to[ito];

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				frame_recv_update(f, isend, jrecv, dx_data, &pat);
				//printf("    dx[%zd, %zd](%zd, %zd) += %g\n",
				//       isend, jrecv, 0, intvl1, dx_data[0]);
			}

			dx_index[0] = coix;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (hsend == coisend || hsend == cojrecv || coisend == cojrecv)
					continue;

				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
				//printf("    dx[%zd, %zd](%zd, %zd) += %g\n",
				//       coisend, cojrecv, intvl1, 0, dx_data[0]);
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
	size_t nintvl1 = nintvl + 1;
	size_t hsend = msg->from;
	size_t dyn_index = v->dyn_index;

	size_t ito, nto = msg->nto;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 };
	struct vpattern pat = vpattern_make(dx_index, 2);
	const struct frame_actor *fa = &f->senders[hsend];
	size_t iz, nz = fa->active.nz;
	const size_t *nmsg;

	//printf("advance { %zd -> [", msg->from);
	//for (ito = 0; ito < nto; ito++) {
	//	printf(" %zd", msg->to[ito]);
	//}
	//printf(" ] } to %zd\n", intvl);


	for (iz = 0, nmsg = fa->nmsg; iz < nz; iz++) {
		size_t jrecv = fa->active.indx[iz];
		size_t coisend = jrecv;
		size_t intvl1;

		for (intvl1 = 0; intvl1 < nintvl1; intvl1++, nmsg++) {
			if (*nmsg == 0)
				continue;

			size_t ix0 = dyn_index + (intvl - 1) + intvl1 * nintvl1;
			size_t ix1 = ix0 + 1;
			size_t coix0 = dyn_index + intvl1 + (intvl - 1) * nintvl1;
			size_t coix1 = coix0 + nintvl1;

			dx_data[0] = -(double)(*nmsg);
			dx_data[1] = +(double)(*nmsg);
			dx_index[0] = ix0;
			dx_index[1] = ix1;

			for (ito = 0; ito < nto; ito++) {
				size_t isend = msg->to[ito];

				if (hsend == isend || hsend == jrecv || isend == jrecv)
					continue;

				frame_recv_update(f, isend, jrecv, dx_data, &pat);
				//printf("    dx[%zd, %zd](%zd, %zd) += %g\n",
				//       isend, jrecv, intvl - 1, intvl1, dx_data[0]);
				//printf("    dx[%zd, %zd](%zd, %zd) += %g\n",
				//       isend, jrecv, intvl, intvl1, dx_data[1]);
			}

			dx_index[0] = coix0;
			dx_index[1] = coix1;

			for (ito = 0; ito < nto; ito++) {
				size_t cojrecv = msg->to[ito];

				if (hsend == coisend || hsend == cojrecv || coisend == cojrecv)
					continue;

				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
				//printf("    dx[%zd, %zd](%zd, %zd) += %g\n",
				//       coisend, cojrecv, intvl1, intvl - 1, dx_data[0]);
				//printf("    dx[%zd, %zd](%zd, %zd) += %g\n",
				//       coisend, cojrecv, intvl1, intvl, dx_data[1]);
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
