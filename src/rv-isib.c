#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *   /----> i
 *  k
 *   \----> j
 */

static char *isib_names[] = { "ISib" };

static void isib_init(struct design_var *v, const struct design *d,
		      void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = isib_names;
}

static void isib_message_add(void *udata, struct frame *f,
			     const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t hsend = msg->from;
	size_t dyn_index = v->dyn_index;
	size_t ito, nto = msg->nto;
	
	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);

	
	const struct frame_actor *fa = &f->senders[hsend];
	size_t iz, nz = fa->active.nz;
	const size_t *nmsg;
	
	for (iz = 0, nmsg = fa->nmsg; iz < nz; iz++) {
		size_t jrecv = fa->active.indx[iz];
		size_t coisend = jrecv;
		
		for (ito = 0; ito < nto; ito++) {
			size_t isend = msg->to[ito];
				
			if (hsend == isend || hsend == jrecv || isend == jrecv)
				continue;
				
			const struct vector *dx = frame_recv_dx(f, isend, jrecv);
			if (vector_item(dx, dyn_index) == 0.0) {
				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
			
		for (ito = 0; ito < nto; ito++) {
			size_t cojrecv = msg->to[ito];
				
			if (hsend == coisend || hsend == cojrecv || coisend == cojrecv)
				continue;
			
			const struct vector *dx = frame_recv_dx(f, coisend, cojrecv);
			if (vector_item(dx, dyn_index) == 0.0) {
				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
			}
		}
	}
}

static struct var_type RECV_VAR_ISIB_REP = {
	VAR_RECV_VAR,
	isib_init,
	NULL,			// deinit
	{
	 isib_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_ISIB = &RECV_VAR_ISIB_REP;
