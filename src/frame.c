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

static void recv_frame_init(struct recv_frame *rf, struct frame *f)
{
	assert(rf);

	intmap_init(&rf->jrecv_dxs, sizeof(struct vector),
		    alignof(struct vector));
	rf->frame = f;
}

static void recv_frame_deinit(struct recv_frame *rf)
{
	assert(rf);

	struct intmap_iter it;
	struct vector *dx;

	INTMAP_FOREACH(it, &rf->jrecv_dxs) {
		dx = INTMAP_VAL(it);
		vector_deinit(dx);
	}
	intmap_deinit(&rf->jrecv_dxs);
}

static void recv_frame_clear(struct recv_frame *rf)
{
	assert(rf);

	struct intmap_iter it;
	struct vector *dx;

	INTMAP_FOREACH(it, &rf->jrecv_dxs) {
		dx = INTMAP_VAL(it);
		vector_deinit(dx);
	}
	intmap_clear(&rf->jrecv_dxs);
}

static struct vector *recv_frame_dx(struct recv_frame *rf, ssize_t jrecv)
{
	assert(rf);

	struct intmap_pos pos;
	struct vector *dx;

	if ((dx = intmap_find(&rf->jrecv_dxs, jrecv, &pos))) {
		return dx;
	}

	dx = intmap_insert(&rf->jrecv_dxs, &pos, NULL);
	vector_init(dx, design_recv_dyn_dim(rf->frame->design));
	return dx;
}

static void var_frames_init(struct frame *frame, const struct design *design)
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

	array_deinit(&frame->vars);
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

static void recv_frames_init(struct frame *f)
{
	assert(f);
	intmap_init(&f->recv_frames, sizeof(struct recv_frame),
		    alignof(struct recv_frame));
}

static void recv_frames_deinit(struct frame *f)
{
	assert(f);
	struct intmap_iter it;
	struct recv_frame *rf;

	INTMAP_FOREACH(it, &f->recv_frames) {
		rf = INTMAP_VAL(it);
		recv_frame_deinit(rf);
	}
	intmap_deinit(&f->recv_frames);
}

static void recv_frames_clear(struct frame *f)
{
	assert(f);
	struct intmap_iter it;
	struct recv_frame *rf;

	INTMAP_FOREACH(it, &f->recv_frames) {
		rf = INTMAP_VAL(it);
		recv_frame_clear(rf);
	}
}

void frame_init(struct frame *f, const struct design *design)
{
	assert(f);
	assert(design);

	f->design = design;
	history_init(&f->history);
	recv_frames_init(f);
	var_frames_init(f, design);
	array_init(&f->events, sizeof(struct frame_event));
	pqueue_init(&f->future_events, frame_event_rcompare,
		    sizeof(struct frame_event));
	array_init(&f->observers, sizeof(struct frame_observer));
	matrix_init(&f->send_xt, design_send_dim(design),
		    design_send_count(design));
	refcount_init(&f->refcount);

	frame_clear(f);
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
	matrix_deinit(&f->send_xt);
	array_deinit(&f->observers);
	pqueue_deinit(&f->future_events);
	array_deinit(&f->events);
	var_frames_deinit(f);
	recv_frames_deinit(f);
	history_deinit(&f->history);
}

void frame_clear(struct frame *f)
{
	assert(f);

	history_clear(&f->history);
	recv_frames_clear(f);
	var_frames_clear(f);
	array_clear(&f->events);
	pqueue_clear(&f->future_events);
	f->next_event_id = 0;

	const struct actors *senders = design_senders(f->design);
	ssize_t ieff = design_send_effects_index(f->design);
	ssize_t itraits = design_send_traits_index(f->design);
	ssize_t ntraits = actors_dim(senders);
	ssize_t nsend = design_send_count(f->design);
	ssize_t isend;

	matrix_fill(&f->send_xt, 0.0);

	for (isend = 0; isend < nsend; isend++) {
		const struct vector *xi = actors_traits(senders, isend);
		struct vector col = matrix_col(&f->send_xt, isend);
		struct vector dst = vector_slice(&col, itraits, ntraits);

		vector_assign_copy(&dst, xi);

		if (design_send_effects(f->design)) {
			matrix_set_item(&f->send_xt, ieff + isend, isend, 1.0);
		}
	}

	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {
		if (obs->h.handle_clear) {
			obs->h.handle_clear(obs->udata, f);
		}
	}
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

static struct recv_frame *frame_recv_frame(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));

	struct intmap_pos pos;
	struct recv_frame *rf;

	if (!(rf = intmap_find(&f->recv_frames, isend, &pos))) {
		rf = intmap_insert(&f->recv_frames, &pos, NULL);
		recv_frame_init(rf, f);
	}
	return rf;
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

void frame_add_observer(struct frame *f, void *udata,
			const struct frame_handlers *h)
{
	assert(f);
	assert(udata);
	assert(h);

	struct frame_observer *obs = array_add(&f->observers, NULL);
	obs->udata = udata;
	obs->h = *h;
}

static bool observer_equals(const void *val, void *udata)
{
	const struct frame_observer *obs = val;
	return obs->udata == udata;
}

void frame_remove_observer(struct frame *f, void *udata)
{
	assert(f);
	assert(udata);

	ssize_t pos =
	    array_find_last_index(&f->observers, observer_equals, udata);
	if (pos >= 0) {
		array_remove_at(&f->observers, pos);
	}
}

static void notify_observers(struct frame *f, const struct frame_event *e)
{
	assert(f);
	assert(e);
	assert(frame_time(f) == e->time);

	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {
		if (obs->h.event_mask & e->type) {
			obs->h.handle_event(obs->udata, e, f);
		}
	}
}

static void notify_vars(struct frame *f, const struct frame_event *e)
{
	assert(f);
	assert(e);
	assert(frame_time(f) == e->time);

	struct frame_var *v;

	ARRAY_FOREACH(v, &f->vars) {
		if (v->design->type->handle_event
		    && v->design->type->event_mask & e->type) {
			v->design->type->handle_event(v, e, f);
		}
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

static void send_var_event_after(struct frame *f, const struct frame_event *fe)
{
	assert(f);
	assert(fe->time == frame_time(f));
	assert(fe->type == SEND_VAR_EVENT);

	const struct send_var_event_meta *meta = &fe->meta.send_var;
	ssize_t isend = meta->item;
	ssize_t index = meta->index;
	double delta = meta->delta;
	double *ptr = matrix_item_ptr(&f->send_xt, index, isend);
	*ptr += delta;
}

static void recv_var_event_after(struct frame *f, const struct frame_event *fe)
{
	assert(f);
	assert(fe->time == frame_time(f));
	assert(fe->type == RECV_VAR_EVENT);

	const struct design *design = f->design;
	ssize_t offset = design_recv_dyn_index(design);
	const struct recv_var_event_meta *meta = &fe->meta.recv_var;
	double *ptr;
	struct vector *dx;

	dx = frame_recv_dx(f, meta->item.isend, meta->item.jrecv);
	ptr = vector_item_ptr(dx, meta->index - offset);
	*ptr += meta->delta;
}

static void event_before(struct frame *f, const struct frame_event *e)
{
	assert(f);
	assert(e->time == frame_time(f));

	notify_vars(f, e);

	switch (e->type) {
	case DYAD_EVENT_INIT:
	case DYAD_EVENT_MOVE:
		dyad_event_before(f, e);
		break;
	case TRIAD_EVENT_INIT:
	case TRIAD_EVENT_MOVE1:
	case TRIAD_EVENT_MOVE2:
		assert(0 && "Not implemented");
	case SEND_VAR_EVENT:
	case RECV_VAR_EVENT:
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
	case SEND_VAR_EVENT:
		send_var_event_after(f, e);
		break;
	case RECV_VAR_EVENT:
		recv_var_event_after(f, e);
		break;
	}

	notify_observers(f, e);
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

struct vector frame_send_x(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));

	return matrix_col(&f->send_xt, isend);
}

struct vector *frame_recv_dx(const struct frame *f, ssize_t isend,
			     ssize_t jrecv)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(0 <= jrecv && jrecv < design_recv_count(f->design));

	struct recv_frame *rf = frame_recv_frame((struct frame *)f, isend);
	return recv_frame_dx(rf, jrecv);
}

void frame_send_mul(double alpha, enum trans_op trans,
		    const struct frame *f,
		    const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_send_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_send_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_send_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_send_dim(f->design));

	enum trans_op transt =
	    (trans == TRANS_NOTRANS ? TRANS_TRANS : TRANS_NOTRANS);
	matrix_mul(alpha, transt, &f->send_xt, x, beta, y);
}

void frame_send_muls(double alpha, enum trans_op trans,
		     const struct frame *f,
		     const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_send_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_send_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_send_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_send_dim(f->design));

	enum trans_op transt =
	    (trans == TRANS_NOTRANS ? TRANS_TRANS : TRANS_NOTRANS);
	matrix_muls(alpha, transt, &f->send_xt, x, beta, y);
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

	const struct design *design = f->design;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);

	design_recv_mul0(alpha, trans, f->design, isend, x, beta, y);
	if (trans == TRANS_NOTRANS) {
		struct svector ys;

		svector_init(&ys, vector_dim(y));
		struct vector xsub = vector_slice(x, off, dim);
		frame_recv_dmul(alpha, trans, f, isend, &xsub, 0.0, &ys);
		svector_axpy(1.0, &ys, y);
		svector_deinit(&ys);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		struct svector ysubs;

		svector_init(&ysubs, vector_dim(&ysub));
		frame_recv_dmul(alpha, trans, f, isend, x, 0.0, &ysubs);
		svector_axpy(1.0, &ysubs, &ysub);
		svector_deinit(&ysubs);
	}
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

	const struct design *design = f->design;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);

	design_recv_muls0(alpha, trans, f->design, isend, x, beta, y);

	if (trans == TRANS_NOTRANS) {
		struct svector xsub;
		struct svector_iter it;

		svector_init(&xsub, dim);
		SVECTOR_FOREACH(it, x) {
			ssize_t i = SVECTOR_IDX(it) - off;
			if (0 <= i && i < dim) {
				svector_set_item(&xsub, i, SVECTOR_VAL(it));
			}
		}
		frame_recv_dmuls(alpha, trans, f, isend, &xsub, 1.0, y);
		svector_deinit(&xsub);
	} else {
		struct vector ysub = vector_slice(y, off, dim);
		frame_recv_dmuls(alpha, trans, f, isend, x, 1.0, &ysub);
	}
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
	       || vector_dim(x) == design_recv_dyn_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_recv_dyn_dim(f->design));

	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}

	const struct design *design = f->design;
	ssize_t dim = design_recv_dyn_dim(design);
	ssize_t i;

	if (dim == 0)
		return;

	const struct recv_frame *rf = intmap_item(&f->recv_frames, isend);
	if (!rf)
		return;

	struct intmap_iter it;
	ssize_t jrecv;
	const struct vector *dx;

	if (trans == TRANS_NOTRANS) {
		INTMAP_FOREACH(it, &rf->jrecv_dxs) {
			jrecv = INTMAP_KEY(it);
			dx = INTMAP_VAL(it);

			double dot = vector_dot(dx, x);
			double *yjrecv = svector_item_ptr(y, jrecv);
			*yjrecv += alpha * dot;
		}
	} else {
		INTMAP_FOREACH(it, &rf->jrecv_dxs) {
			jrecv = INTMAP_KEY(it);
			dx = INTMAP_VAL(it);

			double xjrecv = *vector_item_ptr(x, jrecv);

			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			for (i = 0; i < dim; i++) {
				double val = jscale * vector_item(dx, i);
				double *ptr = svector_item_ptr(y, i);
				*ptr += val;
			}
			// vector_axpy(jscale, dx, y);
		}
	}
}

void frame_recv_dmuls(double alpha, enum trans_op trans,
		      const struct frame *f, ssize_t isend,
		      const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_recv_dyn_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_recv_count(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_recv_dyn_dim(f->design));

	/* y := beta y */
	if (beta == 0.0) {
		vector_fill(y, 0.0);
	} else if (beta != 1.0) {
		vector_scale(y, beta);
	}

	const struct design *design = f->design;
	ssize_t dim = design_recv_dyn_dim(design);

	if (dim == 0)
		return;

	const struct recv_frame *rf = intmap_item(&f->recv_frames, isend);
	if (!rf)
		return;

	struct intmap_iter it;
	ssize_t jrecv;
	const struct vector *dx;

	if (trans == TRANS_NOTRANS) {
		INTMAP_FOREACH(it, &rf->jrecv_dxs) {
			jrecv = INTMAP_KEY(it);
			dx = INTMAP_VAL(it);

			double dot = svector_dot(x, dx);
			double *yjrecv = vector_item_ptr(y, jrecv);
			*yjrecv += alpha * dot;
		}
	} else {
		INTMAP_FOREACH(it, &rf->jrecv_dxs) {
			jrecv = INTMAP_KEY(it);
			dx = INTMAP_VAL(it);

			double xjrecv = svector_item(x, jrecv);

			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			vector_axpy(jscale, dx, y);
		}
	}
}
