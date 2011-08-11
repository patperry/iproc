#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "ieee754.h"
#include "vars.h"
#include "frame.h"

static void frame_vars_init(struct frame *f)
{
	assert(f);
	assert(frame_design(f));

	const struct design *d = frame_design(f);
	ssize_t i, n = array_count(&d->recv_vars);
	struct design_var *dvs = array_to_ptr(&d->recv_vars);
	struct frame_var *fvs = xcalloc(n, sizeof(*fvs));

	for (i = 0; i < n; i++) {
		struct design_var *dv = &dvs[i];
		struct frame_var *fv = &fvs[i];
		fv->design = dv;

		if (fv->design->type->frame_init) {
			dv->type->frame_init(fv, f);
		}

		frame_add_observer(f, fv, &dv->type->callbacks);
	}
	f->vars = fvs;
}

static void frame_vars_deinit(struct frame *f)
{
	assert(f);

	const struct design *d = frame_design(f);
	ssize_t i, n = array_count(&d->recv_vars);
	struct frame_var *fvs = f->vars;

	for (i = 0; i < n; i++) {
		struct frame_var *fv = &fvs[i];

		frame_remove_observer(f, fv);

		if (fv->design->type->frame_deinit) {
			fv->design->type->frame_deinit(fv);
		}
	}
	xfree(fvs);
}

static void recv_frame_init(struct recv_frame *rf, struct frame *f)
{
	assert(rf);

	array_init(&rf->frame_messages, sizeof(struct frame_message));
	intmap_init(&rf->jrecv_dxs, sizeof(struct vector),
		    alignof(struct vector));
}

static void recv_frame_deinit(struct recv_frame *rf)
{
	assert(rf);

	struct intmap_iter it;
	INTMAP_FOREACH(it, &rf->jrecv_dxs) {
		struct vector *dx = INTMAP_VAL(it);
		vector_deinit(dx);
	}
	intmap_deinit(&rf->jrecv_dxs);

	array_deinit(&rf->frame_messages);
}

static void recv_frame_clear(struct recv_frame *rf)
{
	assert(rf);

	array_clear(&rf->frame_messages);

	struct intmap_iter it;
	INTMAP_FOREACH(it, &rf->jrecv_dxs) {
		struct vector *dx = INTMAP_VAL(it);
		vector_deinit(dx);
	}
	intmap_clear(&rf->jrecv_dxs);
}

static struct vector *recv_frame_dx(struct recv_frame *rf, ssize_t jrecv,
				    ssize_t dyn_dim)
{
	assert(rf);

	struct intmap_pos pos;
	struct vector *dx;

	if (!(dx = intmap_find(&rf->jrecv_dxs, jrecv, &pos))) {
		dx = intmap_insert(&rf->jrecv_dxs, &pos, NULL);
		vector_init(dx, dyn_dim);
	}

	assert(vector_dim(dx) == dyn_dim);
	return dx;
}

static void recv_frames_init(struct frame *f)
{
	assert(f);

	const struct design *d = frame_design(f);
	const struct actors *senders = design_senders(d);
	ssize_t isend, nsend = actors_count(senders);

	struct recv_frame *rfs = xcalloc(nsend, sizeof(struct recv_frame));

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_init(&rfs[isend], f);
	}
	f->recv_frames = rfs;
}

static void recv_frames_deinit(struct frame *f)
{
	assert(f);
	const struct design *d = frame_design(f);
	const struct actors *senders = design_senders(d);
	ssize_t isend, nsend = actors_count(senders);

	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_deinit(&rfs[isend]);
	}
	xfree(rfs);
}

static void recv_frames_clear(struct frame *f)
{
	assert(f);
	const struct design *d = frame_design(f);
	const struct actors *senders = design_senders(d);
	ssize_t isend, nsend = actors_count(senders);

	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_clear(&rfs[isend]);
	}
}

static struct recv_frame *recv_frames_item(const struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(frame_design(f)));

	struct recv_frame *rf = &f->recv_frames[isend];
	return rf;
}

static int frame_event_rcompare(const void *x1, const void *x2)
{
	const struct frame_event *e1 = x1;
	const struct frame_event *e2 = x2;
	return double_rcompare(&e1->time, &e2->time);
}

void frame_init(struct frame *f, const struct design *design)
{
	assert(f);
	assert(design);

	f->design = design;
	f->time = -INFINITY;

	array_init(&f->observers, sizeof(struct frame_observer));
	array_init(&f->current_message_ptrs, sizeof(struct message *));
	pqueue_init(&f->events, frame_event_rcompare,
		    sizeof(struct frame_event));
	recv_frames_init(f);
	frame_vars_init(f);

	frame_clear(f);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	frame_vars_deinit(f);
	recv_frames_deinit(f);
	pqueue_deinit(&f->events);
	array_deinit(&f->current_message_ptrs);
	array_deinit(&f->observers);
}

void frame_clear(struct frame *f)
{
	assert(f);

	f->time = -INFINITY;
	array_clear(&f->current_message_ptrs);
	pqueue_clear(&f->events);
	recv_frames_clear(f);

	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, f);
		}
	}
}

void frame_add_observer(struct frame *f, void *udata,
			const struct frame_callbacks *callbacks)
{
	assert(f);
	assert(udata);
	assert(callbacks);

	struct frame_observer *obs = array_add(&f->observers, NULL);
	obs->udata = udata;
	obs->callbacks = *callbacks;
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

double frame_time(const struct frame *f)
{
	assert(f);
	return f->time;
}

double frame_next_time(const struct frame *f)
{
	assert(f);

	if (pqueue_count(&f->events)) {
		const struct frame_event *e = pqueue_top(&f->events);
		return e->time;
	} else {
		return INFINITY;
	}
}

void frame_add(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	array_add(&f->current_message_ptrs, &msg);
}

static void process_current_messages(struct frame *f)
{
	assert(f);

	if (!array_count(&f->current_message_ptrs))
		return;

	const struct design *d = frame_design(f);
	const struct vector *intvls = design_intervals(d);
	double delta = vector_dim(intvls) ? vector_item(intvls, 0) : INFINITY;
	double tnext = frame_time(f) + delta;

	const struct message **msgp;
	ARRAY_FOREACH(msgp, &f->current_message_ptrs) {
		const struct message *msg = *msgp;
		assert(msg->time == frame_time(f));

		// add the message to the appropriate recv_frame
		struct recv_frame *rf = recv_frames_item(f, msg->from);
		ssize_t imsg = array_count(&rf->frame_messages);
		struct frame_message *fmsg =
		    array_add(&rf->frame_messages, NULL);
		fmsg->msg = msg;
		fmsg->interval = 0;

		// create an advance event (if necessary)
		if (isfinite(tnext)) {
			struct frame_event e;
			e.time = tnext;
			e.isend = msg->from;
			e.imsg = imsg;
			pqueue_push(&f->events, &e);
		}
#ifndef NDEBUG
		ssize_t n = array_count(&f->current_message_ptrs);
		double time = frame_time(f);
#endif

		// notify all observers
		const struct frame_observer *obs;
		ARRAY_FOREACH(obs, &f->observers) {
			if (obs->callbacks.message_add) {
				obs->callbacks.message_add(obs->udata, f, msg);

				// make sure observers don't add messages
				// or advance time
				assert(array_count(&f->current_message_ptrs) ==
				       n);
				assert(frame_time(f) == time);
			}
		}
	}
	array_clear(&f->current_message_ptrs);
}

static void process_event(struct frame *f)
{
	assert(pqueue_count(&f->events));
	assert(!array_count(&f->current_message_ptrs));
	assert(frame_time(f) == frame_next_time(f));

	struct frame_event *e = pqueue_top(&f->events);

	// advance the message interval
	ssize_t isend = e->isend;
	ssize_t imsg = e->imsg;
	struct recv_frame *rf = recv_frames_item(f, isend);
	struct frame_message *fmsg = array_item(&rf->frame_messages, imsg);
	const struct message *msg = fmsg->msg;

	ssize_t intvl0 = fmsg->interval;
	ssize_t intvl = intvl0 + 1;
	fmsg->interval = intvl;

	// update the event queue
	const struct design *d = frame_design(f);
	const struct vector *intvls = design_intervals(d);

	if (intvl < vector_dim(intvls)) {
		double delta = vector_item(intvls, intvl);
		double tnext = msg->time + delta;
		e->time = tnext;
		pqueue_update_top(&f->events);
	} else {
		pqueue_pop(&f->events);
	}

#ifndef NDEBUG
	double time = frame_time(f);
#endif

	// notify all observers
	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {
		if (obs->callbacks.message_advance) {
			obs->callbacks.message_advance(obs->udata, f, msg,
						       intvl);

			// make sure observers don't add messages
			// or advance time
			assert(!array_count(&f->current_message_ptrs));
			assert(frame_time(f) == time);
		}
	}
}

void frame_advance(struct frame *f, double time)
{
	assert(f);
	assert(time >= frame_time(f));

	if (frame_time(f) < time)
		process_current_messages(f);

	while (frame_next_time(f) < time) {
		f->time = frame_next_time(f);
		process_event(f);
	}

	assert(frame_time(f) <= time);
	assert(!array_count(&f->current_message_ptrs));
	assert(frame_next_time(f) >= time);
	f->time = time;
}

void frame_recv_update(struct frame *f, ssize_t isend, ssize_t jrecv,
		       ssize_t dyn_index, double delta)
{
	const struct design *d = frame_design(f);

	assert(0 <= isend && isend < design_send_count(d));
	assert(0 <= jrecv && jrecv < design_recv_count(d));
	assert(0 <= dyn_index && dyn_index < design_recv_dyn_dim(d));
	assert(!isnan(delta));

	struct vector *dx = (struct vector *)frame_recv_dx(f, isend, jrecv);
	double *ptr = vector_item_ptr(dx, dyn_index);
	*ptr += delta;

	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {
		if (obs->callbacks.recv_update) {
			obs->callbacks.recv_update(obs->udata, f, isend,
						   jrecv, dyn_index, delta);
		}
	}
}

//struct vector frame_send_x(struct frame *f, ssize_t isend)
//{
//      assert(f);
//      assert(0 <= isend && isend < design_send_count(f->design));
//
//      return matrix_col(&f->send_xt, isend);
//}

const struct vector *frame_recv_dx(const struct frame *f, ssize_t isend,
				   ssize_t jrecv)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(0 <= jrecv && jrecv < design_recv_count(f->design));

	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);
	const struct design *d = frame_design(f);
	ssize_t dyn_dim = design_recv_dyn_dim(d);
	return recv_frame_dx(rf, jrecv, dyn_dim);
}

/*
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
 */

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

	design_recv_mul0(alpha, trans, f->design, x, beta, y);
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

	design_recv_muls0(alpha, trans, f->design, x, beta, y);

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

	const struct recv_frame *rf = recv_frames_item(f, isend);
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

	const struct recv_frame *rf = recv_frames_item(f, isend);
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
