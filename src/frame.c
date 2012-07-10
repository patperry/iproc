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


static void recv_frame_init(struct recv_frame *rf, struct frame *f)
{
	(void)f;		// unused;
	assert(rf);

	vpattern_init(&rf->active);
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
	vpattern_deinit(&rf->active);
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

static void frame_history_message_add (void *udata, struct history * h, const struct message * msg)
{
	struct frame *f = udata;
	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.message_add) {
			obs->callbacks.message_add(obs->udata, f, msg);
		}
	}
}

static void frame_history_message_advance (void *udata, struct history * h, const struct message * msg, size_t intvl)
{
	struct frame *f = udata;
	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.message_advance) {
			obs->callbacks.message_advance(obs->udata, f, msg, intvl);
		}
	}
}

static void frame_history_clear(void *udata, struct history * h)
{
	
}

static struct history_callbacks frame_history_callbacks = {
	frame_history_message_add,
	frame_history_message_advance,
	frame_history_clear
};


static void frame_dyad_design_update (void *udata, struct design * d, const struct var *v, size_t i, const double *delta,
				      const struct vpattern *pat)
{
	struct frame *f = udata;
	
	struct dyad dyad = frame_ix_dyad(f, i);
	
	int ins;
	vpattern_search(&f->recv_frames[dyad.isend].active, dyad.jrecv, &ins); // add jrecv to active set
	
	/* TODO: pattern is wrong, since it is relative to var */
	
	size_t io, no = f->nobs;
	const struct frame_observer *obs;
	for (io = 0; io < no; io++) {
		obs = &f->observers[io];
		if (obs->callbacks.dyad_update) {
			obs->callbacks.dyad_update(obs->udata, f, dyad.isend, dyad.jrecv, delta, pat);
		}
	}
}


static struct design_callbacks frame_dyad_design_callbacks = {
	frame_dyad_design_update,
	NULL
};


void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl)
{
	assert(f);
	assert(intvls || !nintvl);

	f->nsend = nsend;
	f->nrecv = nrecv;
	f->has_loops = has_loops;
	history_init(&f->history, nsend, nrecv, intvls, nintvl);
	history_add_observer(&f->history, f, &frame_history_callbacks);
	design_init(&f->send_design, f, nsend);
	design_init(&f->recv_design, f, nrecv);
	design_init(&f->dyad_design, f, nsend * nrecv);
	design_add_observer(&f->dyad_design, f, &frame_dyad_design_callbacks);
	
	f->observers = NULL;
	f->nobs = 0;
	f->nobs_max = 0;
	recv_frames_init(f);
	frame_clear(f);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	recv_frames_deinit(f);
	free(f->observers);

	design_remove_observer(&f->dyad_design, f);
	design_deinit(&f->dyad_design);
	
	design_deinit(&f->recv_design);
	design_deinit(&f->send_design);

	history_remove_observer(&f->history, f);
	history_deinit(&f->history);
}

void frame_clear(struct frame *f)
{
	assert(f);

	history_clear(&f->history);
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

void frame_add(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	history_add(&f->history, msg);
}


void frame_advance(struct frame *f, double time)
{
	assert(f);
	history_advance(&f->history, time);
}

const double *frame_recv_dx(const struct frame *f, size_t isend, size_t jrecv)
{
	assert(f);
	assert(isend < frame_send_count(f));
	assert(jrecv < frame_recv_count(f));

	size_t ix = frame_dyad_ix(f, isend, jrecv);
	return design_tvars(&f->dyad_design, ix);
}

void frame_recv_get_dx(const struct frame *f, size_t isend,
		       const double **dxp, const size_t **activep,
		       size_t *nactivep)
{
	assert(isend < frame_send_count(f));

	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);

	*dxp = frame_recv_dx(f, isend, rf->active.indx[0]);
	*activep = rf->active.indx;
	*nactivep = rf->active.nz;
}

