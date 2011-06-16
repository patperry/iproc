#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "ieee754.h"
#include "frame.h"

void fprintf_event(FILE * stream, const struct frame_event *e)
{
	assert(e);

	fprintf(stream, "{\n");
	fprintf(stream, "  id: %" SSIZE_FMT "\n", e->id);
	fprintf(stream, "  time: %" SSIZE_FMT "\n", (ssize_t)(e->time));
	fprintf(stream, "  type: %s\n",
		e->type == DYAD_EVENT_INIT ? "DYAD_EVENT_INIT" :
		e->type == DYAD_EVENT_MOVE ? "DYAD_EVENT_MOVE" :
		e->type == TRIAD_EVENT_INIT ? "TRIAD_EVENT_INIT" :
		e->type == TRIAD_EVENT_MOVE1 ? "TRIAD_EVENT_MOVE1" :
		e->type == TRIAD_EVENT_MOVE2 ? "TRIAD_EVENT_MOVE2" :
		e->type == SEND_VAR_EVENT ? "SEND_VAR_EVENT" :
		e->type == RECV_VAR_EVENT ? "RECV_VAR_EVENT" : "(undefined)");
	fprintf(stream, "}\n");
}

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
	svector_init(dx, design_recv_dim(sf->frame->design));
	return dx;
}

static void var_frames_init(struct frame *frame, struct design *design)
{
	assert(frame);
	assert(design);

	array_init(&frame->vars, sizeof(struct frame_var));
	array_set_capacity(&frame->vars, array_count(&design->recv_vars));

	array_add_range(&frame->vars, NULL, array_count(&design->recv_vars));

	ssize_t i, n = array_count(&design->recv_vars);
	struct design_var *dv;
	struct frame_var *fv;

	for (i = 0; i < n; i++) {
		dv = array_item(&design->recv_vars, i);
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
	array_init(&f->events, sizeof(struct frame_event));
	pqueue_init(&f->future_events, frame_event_rcompare,
		    sizeof(struct frame_event));
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
	pqueue_deinit(&f->future_events);
	array_deinit(&f->events);
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
	array_clear(&f->events);
	pqueue_clear(&f->future_events);
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

static struct send_frame *frame_send_frame(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));

	struct intmap_pos pos;
	struct send_frame *sf;

	if (!(sf = intmap_find(&f->send_frames, isend, &pos))) {
		sf = intmap_insert(&f->send_frames, &pos, NULL);
		send_frame_init(sf, f);
	}
	return sf;
}

ssize_t frame_events_count(const struct frame *f)
{
	assert(f);
	return array_count(&f->events);
}

struct frame_event *frame_events_item(const struct frame *f, ssize_t i)
{
	assert(f);
	assert(0 <= i && i < frame_events_count(f));
	return array_item(&f->events, i);
}

static void notify_listeners(struct frame *f, const struct frame_event *e)
{
	assert(f);
	assert(e);
	assert(frame_time(f) == e->time);

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
}

static void dyad_event_before(struct frame *f, const struct frame_event *fe)
{
	assert(f);
	assert(fe->time == frame_time(f));
	assert(fe->type & (DYAD_EVENT_INIT | DYAD_EVENT_MOVE));

	const struct dyad_event_meta *meta = &fe->meta.dyad;
	double t0, dt;

	if (meta->intvl == vector_dim(&f->design->intervals)) {
		/* pass */
	} else {
		struct frame_event e = *fe;
		t0 = meta->msg_time;
		dt = vector_item(&f->design->intervals, meta->intvl);

		e.time = t0 + dt;
		e.type = DYAD_EVENT_MOVE;
		/* preserve e.id */
		e.meta.dyad.intvl++;

		frame_events_add(f, &e);
	}
}

static void recv_var_event_after(struct frame *f,
				 const struct frame_event *fe)
{
	assert(f);
	assert(fe->time == frame_time(f));
	assert(fe->type == RECV_VAR_EVENT);

	const struct dyad_var_event_meta *meta = &fe->meta.dyad_var;
	struct svector_pos pos;
	double *ptr;
	struct svector *dx;

	dx = frame_recv_dx(f, meta->item.isend, meta->item.jrecv);
	if ((ptr = svector_find(dx, meta->index, &pos))) {
		*ptr += meta->delta;
		if (*ptr == 0.0)
			svector_remove_at(dx, &pos);
	} else {
		svector_insert(dx, &pos, meta->delta);
	}
}

static void event_before(struct frame *f, const struct frame_event *e)
{
	assert(f);
	assert(e->time == frame_time(f));

	notify_listeners(f, e);

	switch (e->type) {
	case DYAD_EVENT_INIT:
	case DYAD_EVENT_MOVE:
		dyad_event_before(f, e);
		break;
	case TRIAD_EVENT_INIT:
	case TRIAD_EVENT_MOVE1:
	case TRIAD_EVENT_MOVE2:
		assert(0 && "Not implemented");
	case RECV_VAR_EVENT:
	case SEND_VAR_EVENT:
		break;
	}
}

static void event_after(struct frame *f, const struct frame_event *e)
{
	assert(f);
	assert(e->time == frame_time(f));
	
	switch (e->type) {
		case DYAD_EVENT_INIT:
		case DYAD_EVENT_MOVE:
		case TRIAD_EVENT_INIT:
		case TRIAD_EVENT_MOVE1:
		case TRIAD_EVENT_MOVE2:
			break;
		case RECV_VAR_EVENT:
			recv_var_event_after(f, e);
			break;
		case SEND_VAR_EVENT:
			assert(0 && "Not implemented");
			break;
	}
}


ssize_t frame_events_add(struct frame *f, struct frame_event *e)
{
	assert(f);
	assert(e);

	if (e->id < 0)
		e->id = f->next_event_id++;

	if (e->time == frame_time(f)) {
		array_add(&f->events, e);
		event_before(f, e);
	} else {
		pqueue_push(&f->future_events, e);
	}

	return e->id;
}

double frame_next_change(const struct frame *f)
{
	assert(f);

	if (pqueue_count(&f->future_events)) {
		const struct frame_event *e = pqueue_top(&f->future_events);
		return e->time;
	} else {
		return INFINITY;
	}
}

void frame_advance(struct frame *f)
{
	assert(f);
	assert(frame_time(f) < INFINITY);

	double t = frame_next_change(f);
	const struct frame_event *e;

	ARRAY_FOREACH(e, &f->events) {
		event_after(f, e);
	}
	array_clear(&f->events);

	history_advance_to(&f->history, t);

	while (pqueue_count(&f->future_events)
	       && (e = pqueue_top(&f->future_events))->time == t) {
		// need to copy from top; process_event may invalidate pe */
		struct frame_event ecopy = *e;
		array_add(&f->events, &ecopy);
		event_before(f, &ecopy);
		pqueue_pop(&f->future_events);
	}
	
	assert(frame_time(f) == t);
}

void frame_advance_to(struct frame *f, double t)
{
	assert(f);
	assert(t >= frame_time(f));
	
	while (frame_next_change(f) <= t) {
		frame_advance(f);
	}
	
	if (frame_time(f) < t) {
		const struct frame_event *e;
		ARRAY_FOREACH(e, &f->events) {
			event_after(f, e);
		}
		array_clear(&f->events);
		history_advance_to(&f->history, t);
	}
	
	assert(frame_time(f) == t);
}

void frame_add(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	ssize_t ito, nto = msg->nto;
	struct frame_event e;
		
	e.type = DYAD_EVENT_INIT;
	e.time = msg->time;
	e.meta.dyad.msg_time = msg->time;
	e.meta.dyad.msg_dyad.isend = msg->from;
	e.meta.dyad.msg_attr = msg->attr;
	e.meta.dyad.intvl = 0;
		
	for (ito = 0; ito < nto; ito++) {
		e.id = -1;
		e.meta.dyad.msg_dyad.jrecv = msg->to[ito];
		frame_events_add(f, &e);
	}
		
	assert(history_tcur(&f->history) == msg->time);
	history_add(&f->history, msg->from, msg->to, msg->nto, msg->attr);
}

struct svector *frame_recv_dx(struct frame *f, ssize_t isend, ssize_t jrecv)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(0 <= jrecv && jrecv < design_recv_count(f->design));

	struct send_frame *sf = frame_send_frame(f, isend);
	return send_frame_dx(sf, jrecv);
}

void frame_recv_mul(double alpha, enum trans_op trans,
	       const struct frame *f, ssize_t isend,
	       const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_recv_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_recv_dim(f->design));

	struct svector diffprod;

	svector_init(&diffprod, vector_dim(y));
	frame_recv_dmul(alpha, trans, f, isend, x, 0.0, &diffprod);
	design_recv_mul0(alpha, trans, f->design, isend, x, beta, y);
	svector_axpy(1.0, &diffprod, y);
	svector_deinit(&diffprod);
}

void frame_recv_muls(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_recv_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_recv_dim(f->design));

	struct svector diffprod;

	svector_init(&diffprod, vector_dim(y));
	frame_recv_dmuls(alpha, trans, f, isend, x, 0.0, &diffprod);
	design_recv_muls0(alpha, trans, f->design, isend, x, beta, y);
	svector_axpy(1.0, &diffprod, y);
	svector_deinit(&diffprod);
}

void frame_recv_dmul(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct vector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_recv_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_recv_dim(f->design));

	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}

	struct design *design = f->design;
	ssize_t ndynamic = design->nrdynamic;

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

void frame_recv_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f, ssize_t isend,
		 const struct svector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_recv_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_recv_dim(f->design));

	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}

	struct design *design = f->design;
	ssize_t ndynamic = design->nrdynamic;

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
