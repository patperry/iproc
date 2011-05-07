#include "port.h"
#include <assert.h>

#include "ieee754.h"
#include "frame.h"
#include "vnrecv.h"


/*
static bool insert(const struct dyad_var *dyad_var, const struct message *msg, struct frame *f, ssize_t index)
{
	assert(dyad_var);
	assert(msg);
	assert(f);
	assert(index >= 0);

	struct vnrecv *v = container_of(dyad_var, struct vnrecv, dyad_var);
	double dt, t, t0 = frame_time(f);
	ssize_t iintvl, nintvl= array_size(&v->intvls);
	ssize_t ito, nto = msg->nto;

	// at t0+, increment interval[0]
	t = double_nextup(t0);
	diff.time = t;
	diff.jrecv = msg->from;
	diff.delta = +1.0;
	diff.index = index;
	for (ito = 0; ito < nto; ito++) {
		diff.isend = msg->to[ito];
		frame_add_dyad_event(f, &diff);
	}
	
	// at (t0 + delta[k])+, move weight from interval[k] to interval[k+1]
	for (iintvl = 0; iintvl < nintvl; iintvl++) {
		dt = *(double *)array_at(&v->intvls, iintvl);
		t = double_nextup(t0 + dt);
		diff.time = t;
		diff.jrecv = msg->from;

		for (ito = 0; ito < nto; ito++) {
			diff.isend = msg->to[ito];
			
			diff.delta = -1.0;
			diff.index = index + iintvl;
			frame_add_dyad_event(f, &diff);
			
			diff.delta = +1.0;
			diff.index = index + iintvl + 1;
			frame_add_dyad_event(f, &diff);
		}
	}
	
	return true;
}
*/

bool vnrecv_init(struct vnrecv *v, const double *intvls, ssize_t n)
{
	assert(v);
	assert(intvls || n == 0);
	assert(n >= 0);
	
	if (array_init(&v->intvls, n, sizeof(double))) {
		array_assign_array(&v->intvls, intvls);
		v->dyad_var.dim = n + 1;
		v->dyad_var.dyad_event_mask = 0;
		v->dyad_var.update_send = NULL;
		v->dyad_var.copy_dx_to = NULL;
		return true;
	}

	return false;
}

void vnrecv_deinit(struct vnrecv *v)
{
	assert(v);
	array_deinit(&v->intvls);
}
