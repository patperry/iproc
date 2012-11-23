#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "history.h"


static void grow_ind_array(struct history_actor *a, size_t nintvl, size_t delta);
static size_t search_ind(struct history_actor *a, size_t nintvl, size_t i);
static size_t find_ind(const struct history_actor *a, size_t i);


void grow_ind_array(struct history_actor *a, size_t nintvl, size_t delta)
{
	size_t nz = a->nz;
	size_t nz1 = nz + delta;

	if (needs_grow(nz1, &a->nzmax)) {
		a->ind = xrealloc(a->ind, a->nzmax * sizeof(a->ind[0]));
		a->nmsg = xrealloc(a->nmsg, a->nzmax * nintvl * sizeof(a->nmsg[0]));
	}
}


size_t search_ind(struct history_actor *a, size_t nintvl, size_t i)
{
	const size_t *base = a->ind, *ptr;
	size_t nz;

	for (nz = a->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - a->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	size_t iz = base - a->ind;
	size_t ntail = a->nz - iz;

	grow_ind_array(a, nintvl, 1);
	memmove(a->ind + iz + 1, a->ind + iz, ntail * sizeof(a->ind[0]));
	memmove(a->nmsg + (iz + 1) * nintvl, a->nmsg + iz * nintvl,
		ntail * nintvl * sizeof(a->nmsg[0]));

	a->ind[iz] = i;
	memset(a->nmsg + iz * nintvl, 0, nintvl * sizeof(a->nmsg[0]));
	a->nz++;

	return iz;
}

size_t find_ind(const struct history_actor *a, size_t i)
{
	const size_t *base = a->ind, *ptr;
	size_t nz;

	for (nz = a->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - a->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	return a->nz;
}


static void history_actor_init(struct history_actor *actr)
{
	actr->message_ixs = NULL;
	actr->nix = 0;
	actr->nix_max = 0;
	actr->nmsg = NULL;
	actr->ind = NULL;
	actr->nz = 0;
	actr->nzmax = 0;
}

static void history_actor_clear(struct history_actor *actr)
{
	actr->nix = 0;
	actr->nz = 0;
}

static void history_actor_grow_ixs(struct history_actor *actr, size_t delta)
{
	if (needs_grow(actr->nix + delta, &actr->nix_max)) {
		actr->message_ixs = xrealloc(actr->message_ixs, actr->nix_max * sizeof(actr->message_ixs[0]));
	}
}

static void history_actor_deinit(struct history_actor *actr)
{
	assert(actr);
	free(actr->nmsg);
	free(actr->ind);
	free(actr->message_ixs);
}


static void history_actors_init(struct history_actor **actrsp, size_t n)
{
	assert(actrsp);

	struct history_actor *actrs = xcalloc(n, sizeof(actrs[0]));
	size_t i;
	for (i = 0; i < n; i++) {
		history_actor_init(&actrs[i]);
	}

	*actrsp = actrs;
}

static void history_actors_clear(struct history_actor *actrs, size_t n)
{
	size_t i;
	for (i = 0; i < n; i++) {
		history_actor_clear(&actrs[i]);
	}
}

static void history_actors_deinit(struct history_actor *actrs, size_t n)
{
	size_t i;
	for (i = 0; i < n; i++) {
		history_actor_deinit(&actrs[i]);
	}

	free(actrs);
}

static int history_event_rcompare(const struct pqueue *q, const void *x1,
				const void *x2)
{
	(void)q;
	const struct history_event *e1 = x1;
	const struct history_event *e2 = x2;
	return double_rcompare(&e1->time, &e2->time);
}

void history_init(struct history *h, size_t nsend, size_t nrecv, const double *intvls, size_t nintvl)
{
	assert(h);
	assert(intvls || !nintvl);
#ifndef NDEBUG
	{
		if (nintvl)
			assert(0 < intvls[0]);

		size_t i;
		for (i = 1; i < nintvl; i++) {
			assert(intvls[i - 1] < intvls[i]);
		}
	}
#endif

	h->nsend = nsend;
	h->nrecv = nrecv;
	h->intvls = xmemdup(intvls, nintvl * sizeof(intvls[0]));
	h->nintvl = nintvl;

	h->observers = NULL;
	h->nobs = 0;
	h->nobs_max = 0;

	h->msgs = NULL;
	h->nmsg = 0;
	h->nmsg_max = 0;

	history_actors_init(&h->senders, nsend);
	history_actors_init(&h->receivers, nrecv);
	h->cur_msg = 0;
	pqueue_init(&h->events, sizeof(struct history_event),
		    history_event_rcompare);
	history_clear(h);
}

void history_deinit(struct history *h)
{
	assert(h);

	size_t nsend = history_send_count(h);
	size_t nrecv = history_recv_count(h);

	pqueue_deinit(&h->events);
	history_actors_deinit(h->receivers, nrecv);
	history_actors_deinit(h->senders, nsend);
	free(h->msgs);
	free(h->observers);
	free(h->intvls);
}

void history_clear(struct history *h)
{
	assert(h);

	size_t nsend = history_send_count(h);
	size_t nrecv = history_recv_count(h);

	h->time = -INFINITY;
	h->nmsg = 0;
	history_actors_clear(h->senders, nsend);
	history_actors_clear(h->receivers, nrecv);
	h->cur_msg = 0;
	pqueue_clear(&h->events);

	size_t i, n = h->nobs;
	const struct history_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &h->observers[i];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, h);
		}
	}
}

static void history_observers_grow(struct history *h, size_t delta)
{
	if (needs_grow(h->nobs + delta, &h->nobs_max)) {
		h->observers = xrealloc(h->observers, h->nobs_max * sizeof(h->observers[0]));
	}
}

void history_add_observer(struct history *h, void *udata,
			const struct history_callbacks *callbacks)
{
	assert(h);
	assert(udata);
	assert(callbacks);

	history_observers_grow(h, 1);
	struct history_observer *obs = &h->observers[h->nobs++];
	obs->udata = udata;
	obs->callbacks = *callbacks;
}

void history_remove_observer(struct history *h, void *udata)
{
	assert(h);
	assert(udata);

	size_t i, n = h->nobs;
	for (i = n; i > 0; i--) {
		if (h->observers[i-1].udata == udata) {
			memmove(h->observers + i - 1, h->observers + i,
				(n - i) * sizeof(h->observers[0]));
			h->nobs = n - 1;
			break;
		}
	}
}

static void history_messages_grow(struct history *h, size_t delta)
{
	if (needs_grow(h->nmsg + delta, &h->nmsg_max)) {
		h->msgs = xrealloc(h->msgs, h->nmsg_max * sizeof(h->msgs[0]));
	}
}

void history_add(struct history *h, const struct message *msg)
{
	assert(h);
	assert(msg);
	assert(msg->time == history_time(h));

	//fprintf(stderr, "================= Add (%p) =================\n", msg);

	history_messages_grow(h, 1);
	struct history_message *hmsg = &h->msgs[h->nmsg++];
	hmsg->message = msg;
	hmsg->interval = 0;
}

static int has_pending_messages(const struct history *h)
{
	assert(h->cur_msg <= h->nmsg);
	return h->cur_msg < h->nmsg;
}

static void process_current_messages(struct history *h)
{
	assert(h);

	if (!has_pending_messages(h))
		return;

	const double *intvls = history_intervals(h);
	size_t nintvl = history_interval_count(h);
	double delta = nintvl ? intvls[0] : INFINITY;
	double tnext = history_time(h) + delta;

	const struct history_message *hmsgs = h->msgs;
	size_t imsg, nmsg = h->nmsg;
	for (imsg = h->cur_msg; imsg < nmsg; imsg++) {
		const struct history_message *hmsg = &hmsgs[imsg];
		const struct message *msg = hmsg->message;
		assert(msg->time == history_time(h));

		// add the message to the sender history
		struct history_actor *actr;
		size_t ito, nto = msg->nto;

		actr = &h->senders[msg->from];
		history_actor_grow_ixs(actr, 1);
		actr->message_ixs[actr->nix++] = imsg;

		for (ito = 0; ito < nto; ito++) {
			size_t ix = search_ind(actr, nintvl, msg->to[ito]);
			size_t *ptr = actr->nmsg + ix * nintvl;
			if (nintvl)
				ptr[0]++;
		}

		// add the message to the receiver histories
		for (ito = 0; ito < nto; ito++) {
			actr = &h->receivers[msg->to[ito]];
			history_actor_grow_ixs(actr, 1);
			actr->message_ixs[actr->nix++] = imsg;
			size_t ix = search_ind(actr, nintvl, msg->from);
			size_t *ptr = actr->nmsg + ix * nintvl;
			if (nintvl)
				ptr[0]++;
		}

		// create an advance event (if necessary)
		if (isfinite(tnext)) {
			struct history_event e;
			e.time = tnext;
			e.imsg = imsg;
			pqueue_push(&h->events, &e);
		}
#ifndef NDEBUG
		size_t nmsg1 = h->nmsg;
		double time = history_time(h);
#endif

		//fprintf(stderr, "-> message_add %p, { %zu -> [ ", msg, msg->from);

		//if (msg->nto > 0)
		//      fprintf(stderr, "%zu", msg->to[0]);
		//for (ito = 1; ito < msg->nto; ito++) {
		//      fprintf(stderr, ", %zu", msg->to[ito]);
		//}
		//fprintf(stderr, " ] }\n");

		// notify all observers
		size_t i, n = h->nobs;
		const struct history_observer *obs;
		for (i = 0; i < n; i++) {
			obs = &h->observers[i];
			if (obs->callbacks.message_add) {
				obs->callbacks.message_add(obs->udata, h, msg);

				// make sure observers don't add messages
				// or advance time
				assert(nmsg == nmsg1);
				assert(history_time(h) == time);
			}
		}
	}
	h->cur_msg = nmsg;
}

static void process_event(struct history *h)
{
	assert(pqueue_count(&h->events));
	assert(!has_pending_messages(h));
	assert(history_time(h) == history_next_time(h));

	struct history_event *e = pqueue_top(&h->events);

	// advance the message interval
	size_t imsg = e->imsg;
	struct history_message *hmsg = history_messages_item(h, imsg);
	const struct message *msg = hmsg->message;

	size_t intvl0 = hmsg->interval;
	size_t intvl = intvl0 + 1;
	hmsg->interval = intvl;

	// update the event queue
	const double *intvls = history_intervals(h);
	size_t nintvl = history_interval_count(h);
	double delta = (intvl < nintvl) ? intvls[intvl] : INFINITY;
	double tnext = msg->time + delta;

	if (isfinite(tnext)) {
		e->time = tnext;
		pqueue_update_top(&h->events);
	} else {
		pqueue_pop(&h->events);
	}

	// update the actors;
	struct history_actor *actr;
	size_t ito, nto = msg->nto;

	actr = &h->senders[msg->from];
	for (ito = 0; ito < nto; ito++) {
		size_t ix = find_ind(actr, msg->to[ito]);
		assert(ix < actr->nz);
		size_t *ptr = actr->nmsg + intvl0 + ix * nintvl;
		
		assert(ptr[0]);
		ptr[0]--;

		if (intvl < nintvl)
			ptr[1]++;
	}

	for (ito = 0; ito < nto; ito++) {
		actr = &h->receivers[msg->to[ito]];
		size_t ix = find_ind(actr, msg->from);
		assert(ix < actr->nz);
		size_t *ptr = actr->nmsg + intvl0 + ix * nintvl;
		
		assert(ptr[0]);
		ptr[0]--;

		if (intvl < nintvl)
			ptr[1]++;
	}

#ifndef NDEBUG
	double time = history_time(h);
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
	size_t i, n = h->nobs;
	const struct history_observer *obs;
	
	for (i = 0; i < n; i++) {
		obs = &h->observers[i];

		if (obs->callbacks.message_advance) {
			obs->callbacks.message_advance(obs->udata, h, msg,
						       intvl);

			// make sure observers don't add messages
			// or advance time
			assert(!has_pending_messages(h));
			assert(history_time(h) == time);
		}
	}
}

void history_advance(struct history *h, double time)
{
	assert(h);
	assert(time >= history_time(h));

	if (history_time(h) < time)
		process_current_messages(h);

	while (history_next_time(h) < time) {
		h->time = history_next_time(h);
		process_event(h);
	}

	assert(history_time(h) <= time);
	assert(!has_pending_messages(h));
	assert(history_next_time(h) >= time);
	h->time = time;
}

