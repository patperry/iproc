#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "design.h"
#include "util.h"
#include "vrecv.h"
#include "frame.h"

uint32_t vrecv_active_hash (const void *x)
{
	const struct vrecv_active *active = x;
	return memory_hash(&active->dyad, sizeof(active->dyad));
}

bool vrecv_active_equals (const void *x, const void *y)
{
	const struct vrecv_active *a = x, *b = y;
	return a->dyad.isend == b->dyad.isend && a->dyad.jrecv == b->dyad.jrecv;
}

static bool vrecv_handle_dyad_event (struct dyad_var *v,
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
	

	struct vrecv *vrecv = container_of(v, struct vrecv, dyad_var);
	const struct vrecv_active *key = container_of(&e->dyad, const struct vrecv_active, dyad);
	struct vrecv_active *active;
	struct svector *dx = frame_dx(f, e->dyad.jrecv, e->dyad.isend);
	
	if (!dx)
		return false;
	
	if (e->type == DYAD_EVENT_INIT) {
		struct hashset_pos pos;
		
		if ((active = hashset_find(&vrecv->active, key, &pos)))
			goto move;
		
		if ((active = hashset_insert(&vrecv->active, &pos, key)))
			goto init;

		return false;
	} else { // e->type == DYAD_EVENT_MOVE
		active = hashset_lookup(&vrecv->active, key);
		assert(active);

		if (active->id == e->id) {
			assert(active->intvl == e->intvl - 1);
			goto move;
		}
		
		return true;
	}

move:		
	svector_clear_at(dx, index + active->intvl);
init:
	if (!svector_set(dx, index + e->intvl, 1.0))
		return false;
	active->id = e->id;
	active->intvl = e->intvl;
	return true;
}


bool vrecv_init(struct vrecv *v, const struct design *d)
{
	assert(v);
	assert(d);
	
	if (!hashset_init(&v->active, vrecv_active_hash, vrecv_active_equals,
			  sizeof(struct vrecv_active)))
	    goto fail_active;
	
	ssize_t n = vector_dim(design_intervals(d));
	v->dyad_var.dim = n + 1;
	v->dyad_var.dyad_event_mask = DYAD_EVENT_INIT | DYAD_EVENT_MOVE;
	v->dyad_var.handle_dyad_event = vrecv_handle_dyad_event;
	return true;
	    
fail_active:
	return false;
}

void vrecv_deinit(struct vrecv *v)
{
	assert(v);
	hashset_deinit(&v->active);
}





static void iproc_vrecip_free(iproc_vrecip * v)
{
	if (v) {
		darray_deinit(&v->intvls);
		free(v);
	}
}

static void design_var_free(iproc_design_var * var)
{
	iproc_vrecip *v = container_of(var, iproc_vrecip, var);
	iproc_vrecip_unref(v);
}

static void
design_var_get_dxs(iproc_design_var * var,
		   iproc_design_ctx * ctx, ssize_t offset)
{
	iproc_vrecip *v = container_of(var, iproc_vrecip, var);
	struct darray *intvls = &v->intvls;
	ssize_t nintvl = darray_size(intvls);
	struct history *history = ctx->history;
	ssize_t isend = ctx->isend;

	if (!history || nintvl == 0)
		return;

	double tcur = history_tcur(history);
	struct event_trace *trace = history_recv(history, isend);

	ssize_t i, n = event_trace_size(trace);
	for (i = 0; i < n; i++) {
		struct events *events = event_trace_at(trace, i);
		
		if (events_empty(events))
			continue;
		
		struct event_meta *meta = events_back(events);
		ssize_t jsend = events_id(events);
		double t = meta->time;
		double dt = tcur - t;
		ssize_t pos = darray_binary_search(intvls,
						   &dt, double_compare);

		if (pos < 0)
			pos = ~pos;

		if (pos < nintvl) {
			/* (jsend, [(pos, +1.0)]) */
			struct svector *dx =
			    iproc_design_ctx_dx(ctx, jsend, false);
			*svector_at(dx, offset + pos) += 1.0;
		}
	}
}

iproc_vrecip *iproc_vrecip_new(double *intvls, ssize_t n)
{
	assert(n >= 0);
	assert(n == 0 || intvls);

	iproc_vrecip *v = calloc(1, sizeof(*v));

	if (!v)
		return NULL;

	iproc_design_var_init(&v->var, n, design_var_get_dxs, design_var_free);
	refcount_init(&v->refcount);

	if (!darray_init(&v->intvls, sizeof(double))) {
		iproc_vrecip_free(v);
		v = NULL;
	} else {
		ssize_t i;
		for (i = 0; i < n; i++) {
			assert(intvls[i] > 0.0);
			assert(i == 0 || intvls[i] > intvls[i - 1]);

			darray_push_back(&v->intvls, intvls + i);
		}
	}

	return v;
}

iproc_vrecip *iproc_vrecip_ref(iproc_vrecip * v)
{
	if (v) {
		refcount_get(&v->refcount);
	}

	return v;
}

static void iproc_vrecip_release(struct refcount *refcount)
{
	iproc_vrecip *v = container_of(refcount, iproc_vrecip, refcount);
	iproc_vrecip_free(v);
}

void iproc_vrecip_unref(iproc_vrecip * v)
{
	if (v) {
		refcount_put(&v->refcount, iproc_vrecip_release);
	}
}
