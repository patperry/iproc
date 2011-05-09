#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "ieee754.h"
#include "frame.h"

static bool send_frame_init(struct send_frame *sf, struct frame *f)
{
	assert(sf);

	if (!intmap_init(&sf->jrecv_dxs, sizeof(struct svector),
			alignof(struct svector)))
		goto fail_jrecv_dxs;

	sf->frame = f;
	return true;

fail_jrecv_dxs:
	return false;
}


static void send_frame_deinit(struct send_frame *sf)
{
	assert(sf);
	
	struct intmap_iter it;
	struct svector *dx;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);
		svector_deinit(dx);
	}
	intmap_iter_deinit(&sf->jrecv_dxs, &it);
	intmap_deinit(&sf->jrecv_dxs);
}


static void send_frame_clear(struct send_frame *sf)
{
	assert(sf);
	
	struct intmap_iter it;
	struct svector *dx;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);
		svector_deinit(dx);
	}
	intmap_iter_deinit(&sf->jrecv_dxs, &it);
	intmap_clear(&sf->jrecv_dxs);
}


static struct svector *send_frame_dx(struct send_frame *sf, ssize_t jrecv)
{
	assert(sf);
	
	struct intmap_pos pos;
	struct svector *dx;
	
	if ((dx = intmap_find(&sf->jrecv_dxs, jrecv, &pos))) {
		return dx;
	}
	
	if ((dx = intmap_insert(&sf->jrecv_dxs, &pos, NULL))) {
		if (svector_init(dx, design_dim(sf->frame->design))) {
			return dx;
		}
		intmap_erase(&sf->jrecv_dxs, &pos);
	}
	
	return NULL;
}


bool frame_init(struct frame *f, struct design *design)
{
	assert(f);
	assert(design);
	
	if (!history_init(&f->history))
	    goto fail_history;
	
	if (!dyad_queue_init(&f->dyad_queue, &design->intervals))
	    goto fail_dyad_queue;
	
	if (!intmap_init(&f->send_frames, sizeof(struct send_frame),
			 alignof(struct send_frame)))
		goto fail_send_frames;
	
	if (!(darray_init(&f->dyad_var_udata, sizeof(void *))
	      && darray_reserve(&f->dyad_var_udata,
				darray_size(&design->design_dyad_vars))))
		goto fail_dyad_var_udata;
	
	darray_resize(&f->dyad_var_udata,
		      darray_size(&design->design_dyad_vars));
	
	ssize_t i, n = darray_size(&design->design_dyad_vars);
	struct design_dyad_var *var;
	void *udata;
	for (i = 0; i < n; i++) {
		var = darray_at(&design->design_dyad_vars, i);
		if (var->var->frame_alloc) {
			if (!(udata = var->var->frame_alloc(var->var, f)))
				goto fail_var_frame_alloc;
			darray_set(&f->dyad_var_udata, i, &udata);
		}
	}
	
	if (!(f->design = design_ref(design)))
		goto fail_design;
	
	if (!refcount_init(&f->refcount))
		goto fail_refcount;
	
	return true;

	refcount_deinit(&f->refcount);
fail_refcount:
	design_free(f->design);
fail_design:
fail_var_frame_alloc:
	for (; i > 0; i--) {
		var = darray_at(&design->design_dyad_vars, i - 1);
		if (var->var->frame_free) {
			udata = *(void **)darray_at(&f->dyad_var_udata, i - 1);
			var->var->frame_free(var->var, f, udata);
		}
	}
	darray_deinit(&f->dyad_var_udata);
fail_dyad_var_udata:
	intmap_deinit(&f->send_frames);
fail_send_frames:
	dyad_queue_deinit(&f->dyad_queue);
fail_dyad_queue:
	history_deinit(&f->history);
fail_history:
	return false;
}


struct frame *frame_alloc(struct design *design)
{
	assert(design);

	struct frame *f;
	
	if (!(f = malloc(sizeof(*f))))
		goto fail_malloc;

	if (!frame_init(f, design))
		goto fail_init;
	
	return f;
	
fail_init:
	free(f);
fail_malloc:
	return NULL;
}

struct frame *frame_ref(struct frame *f)
{
	assert(f);
	
	if (refcount_get(&f->refcount))
		return f;
		
	return NULL;
}

void frame_free(struct frame *f)
{
	if (f && refcount_put(&f->refcount, NULL)) {
		refcount_get(&f->refcount);
		frame_deinit(f);
		free(f);
	}
}

void frame_deinit(struct frame *f)
{
	assert(f);
	
	struct intmap_iter it;
	struct send_frame *sf;

	refcount_deinit(&f->refcount);
	design_free(f->design);
	
	ssize_t i, n = darray_size(&f->design->design_dyad_vars);
	struct design_dyad_var *var;
	void *udata;
	for (i = 0; i < n; i++) {
		var = darray_at(&f->design->design_dyad_vars, i);
		if (var->var->frame_free) {
			udata = *(void **)darray_at(&f->dyad_var_udata, i);
			var->var->frame_free(var->var, f, udata);
		}
	}
	
	darray_deinit(&f->dyad_var_udata);
	
	intmap_iter_init(&f->send_frames, &it);
	while (intmap_iter_advance(&f->send_frames, &it)) {
		sf = intmap_iter_current(&f->send_frames, &it);
		send_frame_deinit(sf);
	}
	intmap_iter_deinit(&f->send_frames, &it);
	intmap_deinit(&f->send_frames);

	dyad_queue_deinit(&f->dyad_queue);
	history_deinit(&f->history);
}

void frame_clear(struct frame *f)
{
	assert(f);
	
	struct intmap_iter it;
	struct send_frame *sf;
	
	ssize_t i, n = darray_size(&f->design->design_dyad_vars);
	struct design_dyad_var *var;
	void *udata;

	for (i = 0; i < n; i++) {
		var = darray_at(&f->design->design_dyad_vars, i);
		if (var->var->frame_clear) {
			udata = *(void **)darray_at(&f->dyad_var_udata, i);
			var->var->frame_clear(var->var, f, udata);
		}
	}
	
	intmap_iter_init(&f->send_frames, &it);
	while (intmap_iter_advance(&f->send_frames, &it)) {
		sf = intmap_iter_current(&f->send_frames, &it);
		send_frame_clear(sf);
	}
	intmap_iter_deinit(&f->send_frames, &it);
	
	history_clear(&f->history);
	dyad_queue_clear(&f->dyad_queue);
}

double frame_time(const struct frame *f)
{
	assert(f);
	return history_tcur(&f->history);
}

double frame_next_update(const struct frame *f)
{
	assert(f);
	
	return dyad_queue_next_update(&f->dyad_queue);
}

const struct history *frame_history(const struct frame *f)
{
	assert(f);
	return &f->history;
}

bool frame_insert(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));
	
	ssize_t i, n = darray_size(&f->design->design_dyad_vars);
	const struct design_dyad_var *design_dyad_var;
	
	for (i = 0; i < n; i++) {
		design_dyad_var = darray_at(&f->design->design_dyad_vars, i);
	}
	
	return (history_insert(&f->history, msg->from, msg->to, msg->nto, msg->attr)
		&& dyad_queue_push(&f->dyad_queue, msg));
}


static struct send_frame *frame_send_frame(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_nsender(f->design));
	
	struct intmap_pos pos;
	struct send_frame *sf;

	if ((sf = intmap_find(&f->send_frames, isend, &pos)))
		return sf;
		
	if ((sf = intmap_insert(&f->send_frames, &pos, NULL))) {
		if (send_frame_init(sf, f))
			return sf;
		
		intmap_erase(&f->send_frames, &pos);
	}

	return NULL;
}


struct svector *frame_dx(struct frame *f, ssize_t isend, ssize_t jrecv)
{
	assert(f);
	assert(0 <= isend && isend < design_nsender(f->design));
	assert(0 <= jrecv && jrecv < design_nreceiver(f->design));
	
	struct send_frame *sf = frame_send_frame(f, isend);
	if (sf) {
		return send_frame_dx(sf, jrecv);
	}
	return NULL;
}


bool frame_advance_to(struct frame *f, double t)
{
	assert(f);
	assert(t >= frame_time(f));
	
	if (t == frame_time(f))
		return true;
	
	if (!history_advance_to(&f->history, t))
		return false;
	
	const struct darray *design_dyad_vars = &f->design->design_dyad_vars;
	struct design_dyad_var *v;
	ssize_t i, n = darray_size(design_dyad_vars);
	const struct dyad_event *e;
	void *udata;
	bool ok;

	while (dyad_queue_next_update(&f->dyad_queue) < t) {
		e = dyad_queue_top(&f->dyad_queue);

		for (i = 0; i < n; i++) {
			v = darray_at(design_dyad_vars, i);
			if (!(v->var->handle_dyad_event
			      && v->var->dyad_event_mask & e->type))
				continue;
			udata = *(void **)darray_at(&f->dyad_var_udata, i);
			ok = v->var->handle_dyad_event(v->var, e, f, v->index, udata);
			if (!ok)
				return false;
		}
		
		dyad_queue_pop(&f->dyad_queue);
	}
	
	return true;
}


bool frame_mul(double alpha, enum trans_op trans,
	       const struct frame *f, ssize_t isend,
	       const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_nsender(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_dim(f->design));
	
	
	struct svector diffprod;
	bool ok = false;
	
	if (svector_init(&diffprod, vector_dim(y))) {
		if (frame_dmul(alpha, trans, f, isend, x, 0.0, &diffprod)) {
			design_mul0(alpha, trans, f->design, isend, x, beta, y);
			svector_axpy(1.0, &diffprod, y);
			ok = true;
		}
		svector_deinit(&diffprod);
	}
	
	return ok;
}


bool frame_muls(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_nsender(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_dim(f->design));
	
	
	struct svector diffprod;
	bool ok = false;
	
	if (svector_init(&diffprod, vector_dim(y))) {
		if (frame_dmuls(alpha, trans, f, isend, x, 0.0, &diffprod)) {
			design_muls0(alpha, trans, f->design, isend, x, beta, y);
			svector_axpy(1.0, &diffprod, y);
			ok = true;
		}
		
		svector_deinit(&diffprod);
	}
	
	return ok;
}


bool frame_dmul(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct vector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_nsender(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_dim(f->design));
	
	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}
	
	struct design *design = f->design;
	ssize_t ndynamic = design->ndynamic;
	
	if (ndynamic == 0)
		return true;

	const struct send_frame *sf = intmap_lookup(&f->send_frames, isend);
	if (!sf)
		return true;
	
	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;
	bool ok = true;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		jrecv = intmap_iter_current_key(&sf->jrecv_dxs, &it);
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);
		if (trans == TRANS_NOTRANS) {
			double dot = svector_dot(dx, x);
			double *yjrecv = svector_at(y, jrecv);
			if (!yjrecv) {
				ok = false;
				break;
			}
			
			*yjrecv += alpha * dot;
		} else {
			double xjrecv = *vector_at(x, jrecv);
			
			if (xjrecv == 0.0)
				continue;
			
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			if (!svector_axpys(jscale, dx, y)) {
				ok = false;
				break;
			}
		}
	}

	intmap_iter_deinit(&sf->jrecv_dxs, &it);
	
	return ok;
}


bool frame_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f, ssize_t isend,
		 const struct svector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_nsender(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_dim(f->design));
	
	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}
	
	struct design *design = f->design;
	ssize_t ndynamic = design->ndynamic;
	
	if (ndynamic == 0)
		return true;
	
	const struct send_frame *sf = intmap_lookup(&f->send_frames, isend);
	if (!sf)
		return true;
	
	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;
	bool ok = true;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		jrecv = intmap_iter_current_key(&sf->jrecv_dxs, &it);
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);

		if (trans == TRANS_NOTRANS) {
			double dot = svector_dots(dx, x);
			double *yjrecv = svector_at(y, jrecv);
			if (!yjrecv) {
				ok = false;
				break;
			}
			
			*yjrecv += alpha * dot;
		} else {
			double xjrecv = svector_get(x, jrecv);
			
			if (xjrecv == 0.0)
				continue;
			
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			if (!svector_axpys(jscale, dx, y)) {
				ok = false;
				break;
			}
		}
	}
	
	intmap_iter_deinit(&sf->jrecv_dxs, &it);
	
	return ok;
}

