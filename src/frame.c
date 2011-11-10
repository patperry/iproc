#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "sblas.h"
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "vars.h"
#include "frame.h"

static void frame_actor_init(struct frame_actor *fa)
{
	assert(fa);
	fa->message_ixs = NULL;
	fa->nix = 0;
	fa->nix_max = 0;
	fa->nmsg = NULL;
	vpattern_init(&fa->active);
}

static void frame_actor_clear(struct frame_actor *fa)
{
	assert(fa);
	fa->nix = 0;
	vpattern_clear(&fa->active);
}

static void frame_actor_grow_ixs(struct frame_actor *fa, size_t delta)
{
	size_t nmax = array_grow(fa->nix, fa->nix_max, delta, SIZE_MAX);
	if (nmax > fa->nix_max) {
		fa->message_ixs = xrealloc(fa->message_ixs, nmax * sizeof(fa->message_ixs[0]));
		fa->nix_max = nmax;
	}
}

static void frame_actor_deinit(struct frame_actor *fa)
{
	assert(fa);
	free(fa->nmsg);
	vpattern_deinit(&fa->active);
	free(fa->message_ixs);
}
static size_t frame_actor_search(struct frame_actor *fa, size_t nintvl1, size_t i)
{
	int ins;
	size_t nzmax = fa->active.nzmax;
	size_t ix = vpattern_search(&fa->active, i, &ins);

	if (ins) {
		if (nzmax != fa->active.nzmax) {
			nzmax = fa->active.nzmax;
			fa->nmsg = xrealloc(fa->nmsg, nzmax * nintvl1 * sizeof(fa->nmsg[0]));
		}
		memmove(fa->nmsg + (ix + 1) * nintvl1, fa->nmsg + ix * nintvl1,
			(fa->active.nz - 1 - ix) * nintvl1 * sizeof(fa->nmsg[0]));
		memset(fa->nmsg + ix * nintvl1, 0, nintvl1 * sizeof(fa->nmsg[0]));
	}

	return ix;
}

static void frame_actors_init(struct frame_actor **fasp, size_t n)
{
	assert(fasp);

	struct frame_actor *fas = xcalloc(n, sizeof(fas[0]));
	size_t i;
	for (i = 0; i < n; i++) {
		frame_actor_init(&fas[i]);
	}

	*fasp = fas;
}

static void frame_actors_clear(struct frame_actor *fas, size_t n)
{
	size_t i;
	for (i = 0; i < n; i++) {
		frame_actor_clear(&fas[i]);
	}
}

static void frame_actors_deinit(struct frame_actor *fas, size_t n)
{
	size_t i;
	for (i = 0; i < n; i++) {
		frame_actor_deinit(&fas[i]);
	}

	free(fas);
}

static void recv_frame_init(struct recv_frame *rf, struct frame *f)
{
	(void)f;		// unused;
	assert(rf);

	vpattern_init(&rf->active);
	rf->dx = NULL;
}

static void recv_frame_clear(struct recv_frame *rf)
{
	assert(rf);
	vpattern_clear(&rf->active);
}

static void recv_frame_deinit(struct recv_frame *rf)
{
	assert(rf);

	recv_frame_clear(rf);
	free(rf->dx);
	vpattern_deinit(&rf->active);
}

static double *recv_frame_dx(struct recv_frame *rf, size_t jrecv,
			     size_t dyn_dim)
{
	assert(rf);

	int ins;
	size_t nzmax = rf->active.nzmax;
	size_t ix = vpattern_search(&rf->active, jrecv, &ins);

	if (ins) {
		if (nzmax != rf->active.nzmax) {
			nzmax = rf->active.nzmax;
			rf->dx = xrealloc(rf->dx, nzmax * dyn_dim * sizeof(rf->dx[0]));
		}

		memmove(rf->dx + (ix + 1) * dyn_dim,
			rf->dx + ix * dyn_dim,
			(rf->active.nz - 1 - ix) * dyn_dim * sizeof(rf->dx[0]));
		memset(rf->dx + ix * dyn_dim, 0, dyn_dim * sizeof(rf->dx[0]));
	}

	return rf->dx + ix * dyn_dim;
}

static void recv_frames_init(struct frame *f)
{
	assert(f);

	size_t isend, nsend = frame_send_count(f);

	struct recv_frame *rfs = xcalloc(nsend, sizeof(struct recv_frame));

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_init(&rfs[isend], f);
	}
	f->recv_frames = rfs;
}

static void recv_frames_deinit(struct frame *f)
{
	assert(f);
	size_t isend, nsend = frame_send_count(f);

	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_deinit(&rfs[isend]);
	}
	free(rfs);
}

static void recv_frames_clear(struct frame *f)
{
	assert(f);
	size_t isend, nsend = frame_send_count(f);
	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_clear(&rfs[isend]);
	}
}

static struct recv_frame *recv_frames_item(const struct frame *f, size_t isend)
{
	assert(f);
	assert(isend < frame_send_count(f));

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

void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl)
{
	assert(f);
	assert(intvls || !nintvl);
#ifndef NDEBUG
	{
		size_t i;
		for (i = 1; i < nintvl; i++) {
			assert(intvls[i - 1] < intvls[i]);
		}
	}
#endif

	design_init(&f->send_design, f, nsend);
	design_init(&f->recv_design, f, nrecv);
	f->has_loops = has_loops;
	f->intvls = xmemdup(intvls, nintvl * sizeof(intvls[0]));
	f->nintvl = nintvl;

	f->observers = NULL;
	f->nobs = 0;
	f->nobs_max = 0;
	f->frame_messages = NULL;
	f->nfmsg = 0;
	f->nfmsg_max = 0;
	frame_actors_init(&f->senders, nsend);
	frame_actors_init(&f->receivers, nrecv);
	f->cur_fmsg = 0;
	pqueue_init(&f->events, sizeof(struct frame_event),
		    frame_event_rcompare);
	recv_frames_init(f);
	frame_clear(f);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	recv_frames_deinit(f);
	pqueue_deinit(&f->events);
	frame_actors_deinit(f->receivers, nrecv);
	frame_actors_deinit(f->senders, nsend);
	free(f->frame_messages);
	free(f->observers);
	free(f->intvls);
	design_deinit(&f->recv_design);
	design_deinit(&f->send_design);
}

void frame_clear(struct frame *f)
{
	assert(f);

	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	f->time = -INFINITY;
	f->nfmsg = 0;
	frame_actors_clear(f->senders, nsend);
	frame_actors_clear(f->receivers, nrecv);
	f->cur_fmsg = 0;
	pqueue_clear(&f->events);
	recv_frames_clear(f);

	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, f);
		}
	}
}

static void frame_observers_grow(struct frame *f, size_t delta)
{
	size_t nmax = array_grow(f->nobs, f->nobs_max, delta, SIZE_MAX);
	if (nmax > f->nobs_max) {
		f->observers = xrealloc(f->observers, nmax * sizeof(f->observers[0]));
		f->nobs_max = nmax;
	}
}

void frame_add_observer(struct frame *f, void *udata,
			const struct frame_callbacks *callbacks)
{
	assert(f);
	assert(udata);
	assert(callbacks);

	frame_observers_grow(f, 1);
	struct frame_observer *obs = &f->observers[f->nobs++];
	obs->udata = udata;
	obs->callbacks = *callbacks;
}

void frame_remove_observer(struct frame *f, void *udata)
{
	assert(f);
	assert(udata);

	size_t i, n = f->nobs;
	for (i = n; i > 0; i--) {
		if (f->observers[i-1].udata == udata) {
			memmove(f->observers + i - 1, f->observers + i,
				(n - i) * sizeof(f->observers[0]));
			f->nobs = n - 1;
			break;
		}
	}
}

static void frame_messages_grow(struct frame *f, size_t delta)
{
	size_t nmax = array_grow(f->nfmsg, f->nfmsg_max, delta, SIZE_MAX);
	if (nmax > f->nfmsg_max) {
		f->frame_messages = xrealloc(f->frame_messages, nmax * sizeof(f->frame_messages[0]));
		f->nfmsg_max = nmax;
	}
}

void frame_add(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	//fprintf(stderr, "================= Add (%p) =================\n", msg);

	frame_messages_grow(f, 1);
	struct frame_message *fmsg = &f->frame_messages[f->nfmsg++];
	fmsg->message = msg;
	fmsg->interval = 0;
}

static int has_pending_messages(const struct frame *f)
{
	assert(f->cur_fmsg <= f->nfmsg);
	return f->cur_fmsg < f->nfmsg;
}

static void process_current_messages(struct frame *f)
{
	assert(f);

	if (!has_pending_messages(f))
		return;

	const double *intvls = frame_intervals(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;
	double delta = nintvl ? intvls[0] : INFINITY;
	double tnext = frame_time(f) + delta;

	const struct frame_message *fmsgs = f->frame_messages;
	size_t imsg, nmsg = f->nfmsg;
	for (imsg = f->cur_fmsg; imsg < nmsg; imsg++) {
		const struct frame_message *fmsg = &fmsgs[imsg];
		const struct message *msg = fmsg->message;
		assert(msg->time == frame_time(f));

		// add the message to the sender history
		struct frame_actor *fa;
		size_t ito, nto = msg->nto;

		fa = = &f->senders[msg->from];
		frame_actor_grow_ixs(fa, 1);
		fa->message_ixs[fa->nix++] = imsg;

		for (ito = 0; ito < nto; ito++) {
			size_t ix = frame_actor_search(fa, nintvl1, msg->to[ito]);
			size_t *ptr = fa->nmsg + ix * nintvl1;
			ptr[0]++;
		}

		// add the message to the receiver histories
		for (ito = 0; ito < nto; ito++) {
			fa = &f->receivers[msg->to[ito]];
			frame_actor_grow_ixs(fa, 1);
			fa->message_ixs[fa->nix++] = imsg;
			size_t ix = frame_actor_search(fa, nintvl1, msg->from);
			size_t *ptr = fa->nmsg + ix * nintvl1;
			ptr[0]++;
		}

		// create an advance event (if necessary)
		if (isfinite(tnext)) {
			struct frame_event e;
			e.time = tnext;
			e.imsg = imsg;
			pqueue_push(&f->events, &e);
		}
#ifndef NDEBUG
		size_t nmsg1 = f->nfmsg;
		double time = frame_time(f);
#endif

		//fprintf(stderr, "-> message_add %p, { %zu -> [ ", msg, msg->from);

		//if (msg->nto > 0)
		//      fprintf(stderr, "%zu", msg->to[0]);
		//for (ito = 1; ito < msg->nto; ito++) {
		//      fprintf(stderr, ", %zu", msg->to[ito]);
		//}
		//fprintf(stderr, " ] }\n");

		// notify all observers
		size_t i, n = f->nobs;
		const struct frame_observer *obs;
		for (i = 0; i < n; i++) {
			obs = &f->observers[i];
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
	size_t imsg = e->imsg;
	struct frame_message *fmsg = frame_messages_item(f, imsg);
	const struct message *msg = fmsg->message;

	size_t intvl0 = fmsg->interval;
	size_t intvl = intvl0 + 1;
	fmsg->interval = intvl;

	// update the event queue
	const double *intvls = frame_intervals(f);
	size_t nintvl = frame_interval_count(f);
	size_t nintvl1 = nintvl + 1;

	if (intvl < nintvl) {
		double delta = intvls[intvl];
		double tnext = msg->time + delta;
		e->time = tnext;
		pqueue_update_top(&f->events);
	} else {
		pqueue_pop(&f->events);
	}


	// update the actors;
	struct frame_actor *fa;
	size_t ito, nto = msg->nto;

	fa = &f->senders[msg->from];
	for (ito = 0; ito < nto; ito++) {
		size_t ix = frame_actor_search(fa, nintvl1, msg->to[ito]);
		size_t *ptr = fa->nmsg + intvl0 + ix * nintvl1;
		assert(ptr[0]);
		ptr[0]--;
		ptr[1]++;
	}

	for (ito = 0; ito < nto; ito++) {
		fa = &f->receivers[msg->to[ito]];
		size_t ix = frame_actor_search(fa, nintvl1, msg->from);
		size_t *ptr = fa->nmsg + intvl0 + ix * nintvl1;
		assert(ptr[0]);
		ptr[0]--;
		ptr[1]++;
	}

#ifndef NDEBUG
	double time = frame_time(f);
#endif

	//fprintf(stderr, "-> message_advance %p, { %zu" -> [ ", msg, msg->from);
	//size_t ito;

	//if (msg->nto > 0)
	//      fprintf(stderr, "%zu", msg->to[0]);
	//for (ito = 1; ito < msg->nto; ito++) {
	//      fprintf(stderr, ", %zu", msg->to[ito]);
	//}
	//fprintf(stderr, " ] }\n");

	// notify all observers
	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];

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

void frame_recv_update(struct frame *f, size_t isend, size_t jrecv,
		       const double *delta, const struct vpattern *pat)
{
	assert(isend < frame_send_count(f));
	assert(jrecv < frame_recv_count(f));

	double *dx = (double *)frame_recv_dx(f, isend, jrecv);
	sblas_daxpyi(1.0, delta, pat, dx);

	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.recv_update) {
			obs->callbacks.recv_update(obs->udata, f, isend,
						   jrecv, delta, pat);
		}
	}
}

const double *frame_recv_dx(const struct frame *f, size_t isend, size_t jrecv)
{
	assert(f);
	assert(isend < frame_send_count(f));
	assert(jrecv < frame_recv_count(f));

	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);
	const struct design *d = frame_recv_design(f);
	size_t dyn_dim = design_dvars_dim(d);
	return recv_frame_dx(rf, jrecv, dyn_dim);
}

void frame_recv_get_dx(const struct frame *f, size_t isend,
		       double **dxp, size_t **activep,
		       size_t *nactivep)
{
	assert(isend < frame_send_count(f));

	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);

	*dxp = rf->dx;
	*activep = rf->active.indx;
	*nactivep = rf->active.nz;
}

void frame_recv_mul(double alpha, enum blas_trans trans,
		    const struct frame *f, size_t isend,
		    const double *x, double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t off = design_dvars_index(d);

	design_mul0(alpha, trans, d, x, beta, y);

	if (trans == BLAS_NOTRANS) {
		frame_recv_dmul(alpha, trans, f, isend, x + off, 1.0, y);
	} else {
		frame_recv_dmul(alpha, trans, f, isend, x, 1.0, y + off);
	}
}

void frame_recv_muls(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, const struct vpattern *pat,
		     double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t off = design_dvars_index(d);
	size_t dim = design_dvars_dim(d);

	design_muls0(alpha, trans, d, x, pat, beta, y);

	if (trans == BLAS_NOTRANS) {
		ptrdiff_t jx0 = vpattern_find(pat, off);
		ptrdiff_t jx1 = vpattern_find(pat, off + dim);
		size_t jz0 = (jx0 < 0 ? ~jx0 : jx0);
		size_t jz1 = (jx1 < 0 ? ~jx1 : jx1);
		size_t jz;

		const struct recv_frame *rf = recv_frames_item(f, isend);
		size_t iz, nz = rf->active.nz;

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = rf->active.indx[iz];
			const double *dx = rf->dx + iz * dim;

			double dot = 0;
			for (jz = jz0; jz < jz1; jz++) {
				size_t i = pat->indx[jz] - off;
				dot += x[jz] * dx[i];
			}

			y[jrecv] += alpha * dot;
		}
	} else {
		frame_recv_dmuls(alpha, trans, f, isend, x, pat, 1.0, y + off);
	}
}

void frame_recv_dmul(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t n = design_count(d);
	size_t dim = design_dvars_dim(d);
	size_t ny = (trans == BLAS_NOTRANS ? n : dim);

	/* y := beta y */
	if (beta == 0.0) {
		memset(y, 0, ny * sizeof(y[0]));
	} else if (beta != 1.0) {
		blas_dscal(ny, beta, y, 1);
	}

	if (dim == 0)
		return;

	const struct recv_frame *rf = recv_frames_item(f, isend);
	size_t iz, nz = rf->active.nz;
	size_t jrecv;
	const double *dx;

	if (trans == BLAS_NOTRANS) {
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;

			double dot = blas_ddot(dim, dx, 1, x, 1);
			y[jrecv] += alpha * dot;
		}
	} else {
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;

			if (x[jrecv] == 0.0)
				continue;

			blas_daxpy(dim, alpha * x[jrecv], dx, 1, y, 1);
		}
	}
}

void frame_recv_dmuls(double alpha, enum blas_trans trans,
		      const struct frame *f, size_t isend,
		      const double *x, const struct vpattern *pat,
		      double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t n = design_count(d);
	size_t dim = design_dvars_dim(d);
	size_t ny = (trans == BLAS_NOTRANS ? n : dim);

	/* y := beta y */
	if (beta == 0.0) {
		memset(y, 0, ny * sizeof(y[0]));
	} else if (beta != 1.0) {
		blas_dscal(ny, beta, y, 1);
	}

	if (dim == 0)
		return;

	const struct recv_frame *rf = recv_frames_item(f, isend);
	size_t iz, nz = rf->active.nz;
	size_t jrecv;
	const double *dx;

	if (trans == BLAS_NOTRANS) {

		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;

			double dot = sblas_ddoti(x, pat, dx);
			y[jrecv] += alpha * dot;
		}
	} else {
		size_t jz, mz = pat->nz;
		for (jz = 0; jz < mz; jz++) {
			size_t jrecv = pat->indx[jz];
			ptrdiff_t ix = vpattern_find(&rf->active, jrecv);
			if (ix < 0)
				continue;

			dx = rf->dx + ix * dim;
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * x[jz];
			blas_daxpy(dim, jscale, dx, 1, y, 1);
		}
	}
}
