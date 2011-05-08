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


static bool vnrecv_handle_dyad_event (struct dyad_var *v,
				      const struct dyad_event *e,
				      struct frame *f,
				      ssize_t index)
{
	assert(v);
	assert(v->dim >= 0);
	assert(e);
	assert(f);
	assert(f->design);	
	assert(index >= 0);
	assert(index <= design_dim(f->design) - v->dim);
	
	struct svector *dx = frame_dx(f, e->dyad.jrecv, e->dyad.isend);
	struct svector_pos pos;
	double *val;
	
	if (!dx)
		return false;
	
	if (e->type == DYAD_EVENT_MOVE) {
		val = svector_find(dx, index + e->intvl - 1, &pos);
		assert(val);
		*val -= 1.0;
		if (*val == 0.0)
			svector_erase(dx, &pos);
	}

	val = svector_at(dx, index + e->intvl);
	if (!val)
		return false;

	*val += 1.0;
	return true;
}

/*
static bool vnrecv_get_jrecv_dxs(struct dyad_var *dyad_var,
				 struct frame *f, ssize_t index)
{
	assert(dyad_var);
	assert(f);
	assert(index >= 0);
	
	const struct vector *intvls = f->design->intervals;
	struct history *history = &f->history;
	ssize_t isend = f->isend;
	
	double tcur = history_tcur(history);
	struct event_trace *trace = history_recv(history, isend);
	
	ssize_t i, n = event_trace_size(trace);
	for (i = 0; i < n; i++) {
		struct events *events = event_trace_at(trace, i);
		ssize_t jsend = events_id(events);
		ssize_t k, K = events_size(events);
		const struct event_meta *meta;
		struct svector *dx;
		double *dx_pos;
		ssize_t pos;
		double t, dt;
		
		for (k = 0; k < K; k++) {
			meta = events_at(events, k);
			
			t = meta->time;
			dt = tcur - t;
			pos = array_binary_search(&intvls->array, &dt, double_compare);
		
			if (pos < 0)
				pos = ~pos;
			
			// (jsend, [(pos, +1.0)])
			if (!((dx = DEPRECATED_frame_dx(f, jsend))))
				return false;
			if (!((dx_pos = svector_at(dx, index + pos))))
				return false;
			
			*dx_pos += 1.0;
		}
	}
	
	return true;
}
*/

bool vnrecv_init(struct vnrecv *v, const struct design *d)
{
	assert(v);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	v->dyad_var.dim = n + 1;
	v->dyad_var.dyad_event_mask = DYAD_EVENT_INIT | DYAD_EVENT_MOVE;
	v->dyad_var.handle_dyad_event = vnrecv_handle_dyad_event;
	return true;
}

void vnrecv_deinit(struct vnrecv *v)
{
	assert(v);
}
