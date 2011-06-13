#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "ieee754.h"
#include "frame.h"

static int frame_event_rcompare(const void *x1, const void *x2)
{
	const struct frame_event *e1 = x1;
	const struct frame_event *e2 = x2;
	return double_rcompare(&e1->time, &e2->time);
}

static void send_frame_init(struct send_frame *sf, struct frame *f)
{
	assert(sf);

	intmap_init(&sf->jrecv_dxs, sizeof(struct svector),
		    alignof(struct svector));
	sf->frame = f;
}

static void send_frame_deinit(struct send_frame *sf)
{
	assert(sf);

	struct intmap_iter it;
	struct svector *dx;

	INTMAP_FOREACH(it, &sf->jrecv_dxs) {
		dx = INTMAP_VAL(it);
		svector_deinit(dx);
	}
	intmap_deinit(&sf->jrecv_dxs);
}

static void send_frame_clear(struct send_frame *sf)
{
	assert(sf);

	struct intmap_iter it;
	struct svector *dx;

	INTMAP_FOREACH(it, &sf->jrecv_dxs) {
		dx = INTMAP_VAL(it);
		svector_deinit(dx);
	}
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

	dx = intmap_insert(&sf->jrecv_dxs, &pos, NULL);
	svector_init(dx, design_dim(sf->frame->design));
	return dx;
}

static void var_frames_init(struct frame *frame, struct design *design)
{
	assert(frame);
	assert(design);

	array_init(&frame->vars, sizeof(struct frame_var));
	array_set_capacity(&frame->vars, array_count(&design->vars));

	array_add_range(&frame->vars, NULL, array_count(&design->vars));

	ssize_t i, n = array_count(&design->vars);
	struct design_var *dv;
	struct frame_var *fv;

	for (i = 0; i < n; i++) {
		dv = array_item(&design->vars, i);
		fv = array_item(&frame->vars, i);
		fv->design = dv;

		if (fv->design->type->frame_init) {
			dv->type->frame_init(fv, frame);
		}
	}
}

static void var_frames_deinit(struct frame *frame)
{
	assert(frame);

	ssize_t i, n = array_count(&frame->vars);
	struct frame_var *fv;

	for (i = 0; i < n; i++) {
		fv = array_item(&frame->vars, i);

		if (fv->design->type->frame_deinit) {
			fv->design->type->frame_deinit(fv);
		}
	}
}

static void var_frames_clear(struct frame *frame)
{
	assert(frame);

	ssize_t i, n = array_count(&frame->vars);
	struct frame_var *fv;

	for (i = 0; i < n; i++) {
		fv = array_item(&frame->vars, i);

		if (fv->design->type->frame_clear) {
			fv->design->type->frame_clear(fv);
		}
	}
}

static void send_frames_init(struct frame *f)
{
	assert(f);
	intmap_init(&f->send_frames, sizeof(struct send_frame),
		    alignof(struct send_frame));
}

static void send_frames_deinit(struct frame *f)
{
	assert(f);
	struct intmap_iter it;
	struct send_frame *sf;

	INTMAP_FOREACH(it, &f->send_frames) {
		sf = INTMAP_VAL(it);
		send_frame_deinit(sf);
	}
	intmap_deinit(&f->send_frames);
}

static void send_frames_clear(struct frame *f)
{
	assert(f);
	struct intmap_iter it;
	struct send_frame *sf;

	INTMAP_FOREACH(it, &f->send_frames) {
		sf = INTMAP_VAL(it);
		send_frame_clear(sf);
	}
}

void frame_init(struct frame *f, struct design *design)
{
	assert(f);
	assert(design);

	f->design = design_ref(design);
	history_init(&f->history);
	send_frames_init(f);
	var_frames_init(f, design);
	pqueue_init(&f->events, frame_event_rcompare, sizeof(struct frame_event));
	refcount_init(&f->refcount);
	f->next_event_id = 0;
}

struct frame *frame_alloc(struct design *design)
{
	assert(design);

	struct frame *f = xcalloc(1, sizeof(*f));
	frame_init(f, design);
	return f;
}

struct frame *frame_ref(struct frame *f)
{
	assert(f);
	refcount_get(&f->refcount);
	return f;
}

void frame_free(struct frame *f)
{
	if (f && refcount_put(&f->refcount, NULL)) {
		refcount_get(&f->refcount);
		frame_deinit(f);
		xfree(f);
	}
}

void frame_deinit(struct frame *f)
{
	assert(f);

	refcount_deinit(&f->refcount);
	pqueue_deinit(&f->events);
	var_frames_deinit(f);
	send_frames_deinit(f);
	history_deinit(&f->history);
	design_free(f->design);
}

void frame_clear(struct frame *f)
{
	assert(f);

	history_clear(&f->history);
	send_frames_clear(f);
	var_frames_clear(f);
	pqueue_clear(&f->events);
	f->next_event_id = 0;
}

double frame_time(const struct frame *f)
{
	assert(f);
	return history_tcur(&f->history);
}

const struct history *frame_history(const struct frame *f)
{
	assert(f);
	return &f->history;
}

void frame_insert(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	history_insert(&f->history, msg->from, msg->to, msg->nto, msg->attr);
	
	struct frame_event e;
	e.type = MESSAGE_EVENT;
	e.time = msg->time;
	e.meta.message.from = msg->from;
	e.meta.message.to = msg->to;
	e.meta.message.nto = msg->nto;
	e.meta.message.attr = msg->attr;
	frame_event_push(f, &e);
}

static struct send_frame *frame_send_frame(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_nsender(f->design));

	struct intmap_pos pos;
	struct send_frame *sf;

	if (!(sf = intmap_find(&f->send_frames, isend, &pos))) {
		sf = intmap_insert(&f->send_frames, &pos, NULL);
		send_frame_init(sf, f);
	}
	return sf;
}

struct svector *frame_dx(struct frame *f, ssize_t isend, ssize_t jrecv)
{
	assert(f);
	assert(0 <= isend && isend < design_nsender(f->design));
	assert(0 <= jrecv && jrecv < design_nreceiver(f->design));

	struct send_frame *sf = frame_send_frame(f, isend);
	return send_frame_dx(sf, jrecv);
}

void frame_advance_to(struct frame *f, double t)
{
	assert(f);
	assert(t >= frame_time(f));

	double tnext;

	while ((tnext = frame_event_next(f)) < t) {
		history_advance_to(&f->history, tnext);
		frame_event_pop(f);
	}
	
	history_advance_to(&f->history, t);
}

double frame_event_next(const struct frame *f)
{
	assert(f);
	
	if (pqueue_count(&f->events)) {
		return ((struct frame_event *)pqueue_top(&f->events))->time;
	} else {
		return INFINITY;
	}
}

ssize_t frame_event_push(struct frame *f, struct frame_event *e)
{
	assert(f);
	assert(e);
	
	if (e->id < 0)
		e->id = f->next_event_id++;
	pqueue_push(&f->events, e);
	return e->id;
}

static void pop_message_event(struct frame *f)
{
	assert(f);
	assert(pqueue_count(&f->events));
	assert(((struct frame_event *)pqueue_top(&f->events))->type == MESSAGE_EVENT);
	
	struct frame_event *fe = pqueue_top(&f->events);
	struct message_event_meta *meta = &fe->meta.message;
	
	ssize_t ito, nto = meta->nto;
	struct frame_event e;
	
	e.type = DYAD_EVENT_INIT;
	e.time = fe->time;
	e.meta.dyad.msg_time = fe->time;
	e.meta.dyad.msg_dyad.isend = meta->from;
	e.meta.dyad.msg_attr = meta->attr;
	e.meta.dyad.intvl = 0;
	
	for (ito = 0; ito < nto; ito++) {
		e.id = -1;
		e.meta.dyad.msg_dyad.jrecv = meta->to[ito];
		frame_event_push(f, &e);
	}
	
	pqueue_pop(&f->events);
}

static void pop_dyad_event(struct frame *f)
{
	assert(f);
	assert(pqueue_count(&f->events));
	assert(((struct frame_event *)pqueue_top(&f->events))->type & (DYAD_EVENT_INIT | DYAD_EVENT_MOVE));
	
	struct frame_event *fe = pqueue_top(&f->events);
	struct dyad_event_meta *meta = &fe->meta.dyad;
	double t0, dt;
		
	if (meta->intvl == vector_dim(&f->design->intervals)) {
		pqueue_pop(&f->events);
	} else {
		t0 = meta->msg_time;
		dt = vector_item(&f->design->intervals, meta->intvl);
		
		fe->time = t0 + dt;
		fe->type = DYAD_EVENT_MOVE;
		meta->intvl++;
		
		pqueue_update_top(&f->events);
	}

}

static void pop_dyad_var_event(struct frame *f)
{
	assert(f);
	assert(pqueue_count(&f->events));
	assert(((struct frame_event *)pqueue_top(&f->events))->type == DYAD_VAR_EVENT);

	struct frame_event *fe = pqueue_top(&f->events);
	struct dyad_var_event_meta *meta = &fe->meta.dyad_var;
	struct svector_pos pos;
	double *ptr;
	struct svector *dx;
	
	dx = frame_dx(f, meta->item.isend, meta->item.jrecv);
	if ((ptr = svector_find(dx, meta->index, &pos))) {
		*ptr += meta->delta;
		if (*ptr == 0.0)
			svector_remove_at(dx, &pos);
	} else {
		svector_insert(dx, &pos, meta->delta);
	}
	
	pqueue_pop(&f->events);
}

void frame_event_pop(struct frame *f)
{
	assert(f);

	struct frame_event *e = pqueue_top(&f->events);
	
	
	const struct array *vars = &f->vars;
	struct frame_var *v;
	ssize_t i, n = array_count(vars);
	
	for (i = 0; i < n; i++) {
		v = array_item(vars, i);
		if (!(v->design->type->handle_event
		      && v->design->type->event_mask & e->type))
			continue;
			
		v->design->type->handle_event(v, e, f);
	}
	
	switch(e->type) {
	case MESSAGE_EVENT:
		pop_message_event(f);
		break;
	case DYAD_EVENT_INIT:
	case DYAD_EVENT_MOVE:
		pop_dyad_event(f);
		break;
	case TRIAD_EVENT_INIT:
	case TRIAD_EVENT_MOVE1:			
	case TRIAD_EVENT_MOVE2:
		assert(0 && "Not implemented");
	case DYAD_VAR_EVENT:
		pop_dyad_var_event(f);
		break;
	case SENDER_VAR_EVENT:
		assert(0 && "Not implemented");
		break;
	}

}

void frame_mul(double alpha, enum trans_op trans,
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

	svector_init(&diffprod, vector_dim(y));
	frame_dmul(alpha, trans, f, isend, x, 0.0, &diffprod);
	design_mul0(alpha, trans, f->design, isend, x, beta, y);
	svector_axpy(1.0, &diffprod, y);
	svector_deinit(&diffprod);
}

void frame_muls(double alpha, enum trans_op trans,
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

	svector_init(&diffprod, vector_dim(y));
	frame_dmuls(alpha, trans, f, isend, x, 0.0, &diffprod);
	design_muls0(alpha, trans, f->design, isend, x, beta, y);
	svector_axpy(1.0, &diffprod, y);
	svector_deinit(&diffprod);
}

void frame_dmul(double alpha, enum trans_op trans,
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
		return;

	const struct send_frame *sf = intmap_item(&f->send_frames, isend);
	if (!sf)
		return;

	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;

	INTMAP_FOREACH(it, &sf->jrecv_dxs) {
		jrecv = INTMAP_KEY(it);
		dx = INTMAP_VAL(it);
		if (trans == TRANS_NOTRANS) {
			double dot = svector_dot(dx, x);
			double *yjrecv = svector_item_ptr(y, jrecv);
			*yjrecv += alpha * dot;
		} else {
			double xjrecv = *vector_item_ptr(x, jrecv);

			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			svector_axpys(jscale, dx, y);
		}
	}
}

void frame_dmuls(double alpha, enum trans_op trans,
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
		return;

	const struct send_frame *sf = intmap_item(&f->send_frames, isend);
	if (!sf)
		return;

	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;

	INTMAP_FOREACH(it, &sf->jrecv_dxs) {
		jrecv = INTMAP_KEY(it);
		dx = INTMAP_VAL(it);

		if (trans == TRANS_NOTRANS) {
			double dot = svector_dots(dx, x);
			double *yjrecv = svector_item_ptr(y, jrecv);
			*yjrecv += alpha * dot;
		} else {
			double xjrecv = svector_item(x, jrecv);

			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			svector_axpys(jscale, dx, y);
		}
	}
}
