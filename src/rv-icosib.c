#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  v---- i
 *  h
 *  ^---- j
 */

static char *icosib_names[] = { "ICosib" };

static void icosib_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = icosib_names;
}

static void icosib_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim <=
			design_dvars_dim(frame_recv_design(f)));

	size_t dyn_index = v->dyn_index;

	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);
	size_t ito, nto = msg->nto;

	const struct frame_actor *fa;
	size_t iz, nz;
	const size_t *nmsg;
	
	size_t isend = msg->from;
	size_t cojrecv = msg->from;
	
	for (ito = 0; ito < nto; ito++) {
		size_t hrecv = msg->to[ito];
		fa = &f->receivers[hrecv];
		nz = fa->active.nz;
		
		for (iz = 0, nmsg = fa->nmsg; iz < nz; iz++) {
			size_t jrecv = fa->active.indx[iz];
			size_t coisend = jrecv;
			
			if (hrecv != isend && hrecv != jrecv && isend != jrecv) {
				const double *dx = frame_recv_dx(f, isend, jrecv);
				
				if (dx[dyn_index] == 0.0) {
					frame_recv_update(f, isend, jrecv, dx_data, &pat);
				}
			}
				
			if (!hrecv != coisend && hrecv != cojrecv && coisend != cojrecv) {
				const double *dx = frame_recv_dx(f, coisend, cojrecv);
				
				if (dx[dyn_index] == 0.0) {
					frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
				}
			}
		}
	}
}

static struct var_type RECV_VAR_ICOSIB_REP = {
	VAR_RECV_VAR,
	icosib_init,
	NULL,			// deinit
	{
	 icosib_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_ICOSIB = &RECV_VAR_ICOSIB_REP;
