#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  i ----> h ----> j
 */

static char *isend2_names[] = { "ISend2" };

static void isend2_init(struct design_var *v, const struct design *d,
			void *params)
{
	(void)d;		// unused
	(void)params;		// unused;
	assert(v);
	assert(d);
	assert(!params);

	v->dim = 1;
	v->names = isend2_names;
}

static void isend2_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct design_var *v = udata;

	assert(v);
	assert(f);
	assert(msg);
	assert(v->dyn_index + v->dim
	       <= design_dvars_dim(frame_recv_design(f)));

	size_t dyn_index = v->dyn_index;
	size_t ito, nto = msg->nto;
	
	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { dyn_index };
	struct vpattern pat = vpattern_make(dx_index, 1);
	
	const struct frame_actor *fa;
	size_t iz, nz;
	size_t isend = msg->from;
	
	for (ito = 0; ito < nto; ito++) {
		size_t hsend = msg->to[ito];
		fa = &f->senders[hsend];
		nz = fa->active.nz;
		
		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = fa->active.indx[iz];
				
			if (hsend == isend || hsend == jrecv || isend == jrecv)
				continue;
			
			const double *dx = frame_recv_dx(f, isend, jrecv);
			
			if (dx[dyn_index] == 0.0) {
				frame_recv_update(f, isend, jrecv, dx_data, &pat);
			}
		}
	}
	
	size_t cohrecv = isend;
	fa = &f->receivers[cohrecv];
	nz = fa->active.nz;
	
	for (iz = 0; iz < nz; iz++) {
		size_t coisend = fa->active.indx[iz];
		
		for (ito = 0; ito < nto; ito++) {
			size_t cojrecv = msg->to[ito];
				
			if (cohrecv == coisend || cohrecv == cojrecv || coisend == cojrecv)
				continue;
			
			const double *dx = frame_recv_dx(f, coisend, cojrecv);
			
			if (dx[dyn_index] == 0.0) {
				frame_recv_update(f, coisend, cojrecv, dx_data, &pat);
			}
		}
	}
}

static struct var_type RECV_VAR_ISEND2_REP = {
	VAR_RECV_VAR,
	isend2_init,
	NULL,			// deinit
	{
	 isend2_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// clear
	 }
};

const struct var_type *RECV_VAR_ISEND2 = &RECV_VAR_ISEND2_REP;
