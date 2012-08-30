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


void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl)
{
	assert(f);
	assert(intvls || !nintvl);

	f->nsend = nsend;
	f->nrecv = nrecv;
	f->has_loops = has_loops;
	f->observers = NULL;
	f->nobs = 0;
	f->nobs_max = 0;

	history_init(&f->history, nsend, nrecv, intvls, nintvl);
	history_add_observer(&f->history, f, &frame_history_callbacks);
	design_init(&f->send_design, f, nsend);
	design_init(&f->recv_design, f, nrecv);
	design2_init(&f->dyad_design, f, nsend, nrecv);
	
	frame_clear(f);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	free(f->observers);

	design2_deinit(&f->dyad_design);
	
	design_deinit(&f->recv_design);
	design_deinit(&f->send_design);

	history_remove_observer(&f->history, f);
	history_deinit(&f->history);
}

void frame_clear(struct frame *f)
{
	assert(f);

	history_clear(&f->history);
	
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

