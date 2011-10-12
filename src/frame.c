#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "sblas.h"
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "vars.h"
#include "frame.h"

static void frame_actor_init(struct frame_actor *fa)
{
	assert(fa);
	array_init(&fa->message_ixs, sizeof(ssize_t));
}

static void frame_actor_clear(struct frame_actor *fa)
{
	assert(fa);
	array_clear(&fa->message_ixs);
}

static void frame_actor_deinit(struct frame_actor *fa)
{
	assert(fa);
	array_deinit(&fa->message_ixs);
}

static void frame_actors_init(struct frame_actor **pfas, ssize_t n)
{
	assert(pfas);
	assert(n >= 0);

	struct frame_actor *fas = xcalloc(n, sizeof(fas[0]));
	ssize_t i;
	for (i = 0; i < n; i++) {
		frame_actor_init(&fas[i]);
	}

	*pfas = fas;
}

static void frame_actors_clear(struct frame_actor *fas, ssize_t n)
{
	assert(fas);
	assert(n >= 0);

	ssize_t i;
	for (i = 0; i < n; i++) {
		frame_actor_clear(&fas[i]);
	}
}

static void frame_actors_deinit(struct frame_actor *fas, ssize_t n)
{
	assert(fas);
	assert(n >= 0);

	ssize_t i;
	for (i = 0; i < n; i++) {
		frame_actor_deinit(&fas[i]);
	}

	free(fas);
}

static void frame_senders_init(struct frame *f)
{
	const struct design *d = frame_design(f);
	ssize_t n = design_send_count(d);
	frame_actors_init(&f->senders, n);
}

static void frame_senders_clear(struct frame *f)
{
	const struct design *d = frame_design(f);
	ssize_t n = design_send_count(d);
	frame_actors_clear(f->senders, n);
}

static void frame_senders_deinit(struct frame *f)
{
	const struct design *d = frame_design(f);
	ssize_t n = design_send_count(d);
	frame_actors_deinit(f->senders, n);
}

static struct frame_actor *frame_senders_item(const struct frame *f, ssize_t i)
{
	assert(f);
	assert(0 <= i && i < frame_senders_count(f));
	return &f->senders[i];
}

static void frame_receivers_init(struct frame *f)
{
	const struct design *d = frame_design(f);
	ssize_t n = design_recv_count(d);
	frame_actors_init(&f->receivers, n);
}

static void frame_receivers_clear(struct frame *f)
{
	const struct design *d = frame_design(f);
	ssize_t n = design_recv_count(d);
	frame_actors_clear(f->receivers, n);
}

static void frame_receivers_deinit(struct frame *f)
{
	const struct design *d = frame_design(f);
	ssize_t n = design_recv_count(d);
	frame_actors_deinit(f->receivers, n);
}

static struct frame_actor *frame_receivers_item(const struct frame *f,
						ssize_t i)
{
	assert(f);
	assert(0 <= i && i < frame_receivers_count(f));
	return &f->receivers[i];
}

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
	free(fvs);
}

static void recv_frame_init(struct recv_frame *rf, struct frame *f)
{
	(void)f;		// unused;
	assert(rf);

	rf->active = NULL;
	rf->dx = NULL;
	rf->nactive = 0;
	rf->nactive_max = 0;
}

static void recv_frame_grow_active(struct recv_frame *rf)
{
	size_t nactive_max = rf->nactive_max;

	if (rf->nactive == nactive_max) {
		nactive_max = ARRAY_GROW(nactive_max, SIZE_MAX);
		rf->dx = xrealloc(rf->dx, nactive_max * sizeof(rf->dx[0]));
		rf->active =
		    xrealloc(rf->active, nactive_max * sizeof(rf->active[0]));
		rf->nactive_max = nactive_max;
	}
}

static void recv_frame_clear(struct recv_frame *rf)
{
	assert(rf);

	size_t i, n = rf->nactive;
	for (i = 0; i < n; i++) {
		vector_deinit(&rf->dx[i]);
	}
	rf->nactive = 0;
}

static void recv_frame_deinit(struct recv_frame *rf)
{
	assert(rf);

	recv_frame_clear(rf);
	free(rf->dx);
	free(rf->active);
}

static struct vector *recv_frame_dx(struct recv_frame *rf, ssize_t jrecv,
				    ssize_t dyn_dim)
{
	assert(rf);

	ptrdiff_t ix = sblas_find(rf->nactive, rf->active, jrecv);
	if (ix < 0) {
		recv_frame_grow_active(rf);
		ix = ~ix;
		memmove(rf->active + ix + 1, rf->active + ix,
			(rf->nactive - ix) * sizeof(rf->active[0]));
		memmove(rf->dx + ix + 1, rf->dx + ix,
			(rf->nactive - ix) * sizeof(rf->dx[0]));
		rf->active[ix] = jrecv;
		vector_init(&rf->dx[ix], dyn_dim);
		rf->nactive++;
	}

	struct vector *dx = &rf->dx[ix];
	assert(vector_dim(dx) == dyn_dim);
	return dx;
}

static void recv_frames_init(struct frame *f)
{
	assert(f);

	const struct design *d = frame_design(f);
	ssize_t isend, nsend = design_send_count(d);

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
	ssize_t isend, nsend = design_send_count(d);

	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_deinit(&rfs[isend]);
	}
	free(rfs);
}

static void recv_frames_clear(struct frame *f)
{
	assert(f);
	const struct design *d = frame_design(f);
	ssize_t isend, nsend = design_send_count(d);
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

static int frame_event_rcompare(const struct pqueue *q, const void *x1,
				const void *x2)
{
	(void)q;
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
	array_init(&f->frame_messages, sizeof(struct frame_message));
	frame_senders_init(f);
	frame_receivers_init(f);
	f->cur_fmsg = 0;
	pqueue_init(&f->events, sizeof(struct frame_event),
		    frame_event_rcompare);

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
	frame_receivers_deinit(f);
	frame_senders_deinit(f);
	array_deinit(&f->frame_messages);
	array_deinit(&f->observers);
}

void frame_clear(struct frame *f)
{
	assert(f);

	f->time = -INFINITY;
	array_clear(&f->frame_messages);
	frame_senders_clear(f);
	frame_receivers_clear(f);
	f->cur_fmsg = 0;
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

void frame_add(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	//fprintf(stderr, "================= Add (%p) =================\n", msg);

	struct frame_message *fmsg = array_add(&f->frame_messages, NULL);
	fmsg->message = msg;
	fmsg->interval = 0;
}

static bool has_pending_messages(const struct frame *f)
{
	assert(0 <= f->cur_fmsg
	       && f->cur_fmsg <= array_count(&f->frame_messages));
	return f->cur_fmsg < array_count(&f->frame_messages);
}

static void process_current_messages(struct frame *f)
{
	assert(f);

	if (!has_pending_messages(f))
		return;

	const struct design *d = frame_design(f);
	const struct vector *intvls = design_intervals(d);
	double delta = vector_dim(intvls) ? vector_item(intvls, 0) : INFINITY;
	double tnext = frame_time(f) + delta;

	const struct frame_message *fmsgs = array_to_ptr(&f->frame_messages);
	ssize_t imsg, nmsg = array_count(&f->frame_messages);
	for (imsg = f->cur_fmsg; imsg < nmsg; imsg++) {
		const struct frame_message *fmsg = &fmsgs[imsg];
		const struct message *msg = fmsg->message;
		assert(msg->time == frame_time(f));

		// add the message to the sender history
		struct frame_actor *fa;
		fa = frame_senders_item(f, msg->from);
		array_add(&fa->message_ixs, &imsg);

		// add the message to the receiver histories
		ssize_t ito, nto = msg->nto;
		for (ito = 0; ito < nto; ito++) {
			fa = frame_receivers_item(f, msg->to[ito]);
			array_add(&fa->message_ixs, &imsg);
		}

		// create an advance event (if necessary)
		if (isfinite(tnext)) {
			struct frame_event e;
			e.time = tnext;
			e.imsg = imsg;
			pqueue_push(&f->events, &e);
		}
#ifndef NDEBUG
		ssize_t nmsg1 = array_count(&f->frame_messages);
		double time = frame_time(f);
#endif

		//fprintf(stderr, "-> message_add %p, { %" SSIZE_FMT " -> [ ", msg, msg->from);

		//if (msg->nto > 0)
		//      fprintf(stderr, "%" SSIZE_FMT, msg->to[0]);
		//for (ito = 1; ito < msg->nto; ito++) {
		//      fprintf(stderr, ", %" SSIZE_FMT, msg->to[ito]);
		//}
		//fprintf(stderr, " ] }\n");

		// notify all observers
		const struct frame_observer *obs;
		ARRAY_FOREACH(obs, &f->observers) {
			if (obs->callbacks.message_add) {
				obs->callbacks.message_add(obs->udata, f, msg);

				// make sure observers don't add messages
				// or advance time
				assert(nmsg == nmsg1);
				assert(frame_time(f) == time);
			}
		}
	}
	f->cur_fmsg = nmsg;
}

static void process_event(struct frame *f)
{
	assert(pqueue_count(&f->events));
	assert(!has_pending_messages(f));
	assert(frame_time(f) == frame_next_time(f));

	struct frame_event *e = pqueue_top(&f->events);

	// advance the message interval
	ssize_t imsg = e->imsg;
	struct frame_message *fmsg = frame_messages_item(f, imsg);
	const struct message *msg = fmsg->message;

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

	//fprintf(stderr, "-> message_advance %p, { %" SSIZE_FMT " -> [ ", msg, msg->from);
	//ssize_t ito;

	//if (msg->nto > 0)
	//      fprintf(stderr, "%" SSIZE_FMT, msg->to[0]);
	//for (ito = 1; ito < msg->nto; ito++) {
	//      fprintf(stderr, ", %" SSIZE_FMT, msg->to[ito]);
	//}
	//fprintf(stderr, " ] }\n");

	// notify all observers
	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {

		if (obs->callbacks.message_advance) {
			obs->callbacks.message_advance(obs->udata, f, msg,
						       intvl);

			// make sure observers don't add messages
			// or advance time
			assert(!has_pending_messages(f));
			assert(frame_time(f) == time);
		}
	}
}

void frame_advance(struct frame *f, double time)
{
	assert(f);
	assert(time >= frame_time(f));

	//fprintf(stderr, "=============== Advance (%"SSIZE_FMT") ===============\n", (ssize_t)time);

	if (frame_time(f) < time)
		process_current_messages(f);

	while (frame_next_time(f) < time) {
		f->time = frame_next_time(f);
		process_event(f);
	}

	assert(frame_time(f) <= time);
	assert(!has_pending_messages(f));
	assert(frame_next_time(f) >= time);
	f->time = time;
}

void frame_recv_update(struct frame *f, ssize_t isend, ssize_t jrecv,
		       const struct svector *delta)
{
#ifndef NDEBUG
	const struct design *d = frame_design(f);
#endif
	assert(0 <= isend && isend < design_send_count(d));
	assert(0 <= jrecv && jrecv < design_recv_count(d));
	assert(svector_dim(delta) == design_recv_dyn_dim(d));

	struct vector *dx = (struct vector *)frame_recv_dx(f, isend, jrecv);
	svector_axpy(1.0, delta, dx);

	//fprintf(stderr, "-> recv_update %" SSIZE_FMT ", %" SSIZE_FMT "\n", isend, jrecv);

	const struct frame_observer *obs;
	ARRAY_FOREACH(obs, &f->observers) {
		if (obs->callbacks.recv_update) {

			obs->callbacks.recv_update(obs->udata, f, isend,
						   jrecv, delta);
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

void frame_recv_get_dx(const struct frame *f, ptrdiff_t isend,
		       struct vector **dxp, ptrdiff_t **activep,
		       size_t *nactivep)
{
	assert(f);
	assert(0 <= isend && isend < design_send_count(f->design));
	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);
	*dxp = rf->dx;
	*activep = rf->active;
	*nactivep = rf->nactive;
}

/*
void frame_send_mul(double alpha, enum blas_trans trans,
		    const struct frame *f,
		    const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || vector_dim(x) == design_send_dim(f->design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_send_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(x) == design_send_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_send_dim(f->design));

	enum blas_trans transt =
	    (trans == BLAS_NOTRANS ? TRANS_TRANS : BLAS_NOTRANS);
	matrix_mul(alpha, transt, &f->send_xt, x, beta, y);
}

void frame_send_muls(double alpha, enum blas_trans trans,
		     const struct frame *f,
		     const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || svector_dim(x) == design_send_dim(f->design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_send_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || svector_dim(x) == design_send_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_send_dim(f->design));

	enum blas_trans transt =
	    (trans == BLAS_NOTRANS ? TRANS_TRANS : BLAS_NOTRANS);
	matrix_muls(alpha, transt, &f->send_xt, x, beta, y);
}
 */

void frame_recv_mul(double alpha, enum blas_trans trans,
		    const struct frame *f, ssize_t isend,
		    const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || vector_dim(x) == design_recv_dim(f->design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(x) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_recv_dim(f->design));

	const struct design *design = f->design;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);

	design_recv_mul0(alpha, trans, f->design, x, beta, y);
	if (trans == BLAS_NOTRANS) {
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

void frame_recv_muls(double alpha, enum blas_trans trans,
		     const struct frame *f, ssize_t isend,
		     const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || svector_dim(x) == design_recv_dim(f->design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || svector_dim(x) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(y) == design_recv_dim(f->design));

	const struct design *design = f->design;
	ssize_t off = design_recv_dyn_index(design);
	ssize_t dim = design_recv_dyn_dim(design);

	design_recv_muls0(alpha, trans, f->design, x, beta, y);

	if (trans == BLAS_NOTRANS) {
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

void frame_recv_dmul(double alpha, enum blas_trans trans,
		     const struct frame *f, ssize_t isend,
		     const struct vector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || vector_dim(x) == design_recv_dyn_dim(f->design));
	assert(trans != BLAS_NOTRANS
	       || svector_dim(y) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || vector_dim(x) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
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
	ssize_t jrecv;
	const struct vector *dx;

	if (trans == BLAS_NOTRANS) {
		size_t iz, nz = rf->nactive;
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active[iz];
			dx = &rf->dx[iz];

			double dot = vector_dot(dx, x);
			double *yjrecv = svector_item_ptr(y, jrecv);
			*yjrecv += alpha * dot;
		}
	} else {
		size_t iz, nz = rf->nactive;
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active[iz];
			dx = &rf->dx[iz];

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

void frame_recv_dmuls(double alpha, enum blas_trans trans,
		      const struct frame *f, ssize_t isend,
		      const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(0 <= isend && isend < design_send_count(f->design));
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS
	       || svector_dim(x) == design_recv_dyn_dim(f->design));
	assert(trans != BLAS_NOTRANS
	       || vector_dim(y) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
	       || svector_dim(x) == design_recv_count(f->design));
	assert(trans == BLAS_NOTRANS
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
	ssize_t jrecv;
	const struct vector *dx;

	if (trans == BLAS_NOTRANS) {
		size_t iz, nz = rf->nactive;
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active[iz];
			dx = &rf->dx[iz];

			double dot = svector_dot(x, dx);
			double *yjrecv = vector_item_ptr(y, jrecv);
			*yjrecv += alpha * dot;
		}
	} else {
		size_t iz, nz = rf->nactive;
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active[iz];
			dx = &rf->dx[iz];

			double xjrecv = svector_item(x, jrecv);

			if (xjrecv == 0.0)
				continue;

			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			vector_axpy(jscale, dx, y);
		}
	}
}
