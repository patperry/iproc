#ifndef HISTORY_H
#define HISTORY_H

#include <assert.h>
#include <math.h>
#include "messages.h"
#include "pqueue.h"

struct history {
	size_t nsend, nrecv;
	double *intvls;
	size_t nintvl;

	double time;

	struct history_observer *observers;
	size_t nobs, nobs_max;
	struct history_message *msgs;
	size_t nmsg, nmsg_max;
	struct history_actor *senders;
	struct history_actor *receivers;
	size_t cur_msg;
	struct pqueue events;
};

struct history_callbacks {
	void (*message_add) (void *udata, struct history * h,
			     const struct message * msg);
	void (*message_advance) (void *udata, struct history * h,
				 const struct message * msg, size_t intvl);
	void (*clear) (void *udata, struct history * h);
};

struct history_message {
	const struct message *message;
	size_t interval;
};

struct history_event {
	double time;
	size_t imsg;
};

struct history_actor {
	size_t *message_ixs;
	size_t nix, nix_max;
	size_t *nmsg;
	size_t *ind;
	size_t nz, nzmax;
};

struct history_observer {
	void *udata;
	struct history_callbacks callbacks;
};

/* create/destroy/clear */
void history_init(struct history *h, size_t nsend, size_t nrecv,
		const double *intvls, size_t nintvl);
void history_deinit(struct history *h);
void history_clear(struct history *h);

/* properties */
static inline const double *history_intervals(const struct history *h);
static inline size_t history_interval_count(const struct history *h);


/* time */
static inline double history_time(const struct history *h);	// current time
static inline double history_next_time(const struct history *h);	// next change
void history_advance(struct history *h, double time);	// advance time

/* messages */
static inline size_t history_messages_count(const struct history *h);
static inline struct history_message *history_messages_item(const struct history *h,
							size_t i);
void history_add(struct history *h, const struct message *msg);

/* actors */
static inline size_t history_send_count(const struct history *h);
static inline size_t history_recv_count(const struct history *h);

static inline void history_get_send_messages(const struct history *h, size_t isend,
					     const size_t **imsg, size_t *nmsg);
static inline void history_get_send_active(const struct history *h, size_t isend,
					   const size_t **jrecv, size_t *nrecv);
static inline const size_t *history_send_counts(const struct history *h, size_t isend);

static inline void history_get_recv_messages(const struct history *h, size_t irecv,
					     const size_t **imsg, size_t *nmsg);
static inline void history_get_recv_active(const struct history *h, size_t jrecv,
					   const size_t **isend, size_t *nsend);
static inline const size_t *history_recv_counts(const struct history *h, size_t jrecv);



/* observers */
void history_add_observer(struct history *h, void *udata,
			  const struct history_callbacks *callbacks);
void history_remove_observer(struct history *h, void *udata);


/* inline function definitions */
const double *history_intervals(const struct history *h)
{
	return h->intvls;
}

size_t history_interval_count(const struct history *h)
{
	return h->nintvl;
}

double history_time(const struct history *h)
{
	assert(h);
	return h->time;
}

double history_next_time(const struct history *h)
{
	assert(h);

	if (pqueue_count(&h->events)) {
		const struct history_event *e = pqueue_top(&h->events);
		return e->time;
	} else {
		return INFINITY;
	}
}

size_t history_messages_count(const struct history *h)
{
	assert(h);
	return h->nmsg;
}

struct history_message *history_messages_item(const struct history *h, size_t imsg)
{
	assert(h);
	assert(imsg < history_messages_count(h));
	return &h->msgs[imsg];
}

size_t history_send_count(const struct history *h)
{
	assert(h);
	return h->nsend;
}

size_t history_recv_count(const struct history *h)
{
	assert(h);
	return h->nrecv;
}

void history_get_send_messages(const struct history *h, size_t isend,
			       const size_t **imsg, size_t *nmsg)
{
	assert(h);
	assert(isend < history_send_count(h));
	assert(imsg);
	assert(nmsg);

	const struct history_actor *a = &h->senders[isend];
	*imsg = a->message_ixs;
	*nmsg = a->nix;
}

void history_get_send_active(const struct history *h, size_t isend, const size_t **jrecv,
			     size_t *nrecv)
{
	const struct history_actor *a = &h->senders[isend];
	*jrecv = a->ind;
	*nrecv = a->nz;
}

const size_t *history_send_counts(const struct history *h, size_t isend)
{
	assert(h);
	assert(isend < history_send_count(h));
	const struct history_actor *a = &h->senders[isend];
	return a->nmsg;
}

void history_get_recv_messages(const struct history *h, size_t irecv,
			       const size_t **imsg, size_t *nmsg)
{
	assert(h);
	assert(irecv < history_recv_count(h));
	assert(imsg);
	assert(nmsg);

	const struct history_actor *a = &h->receivers[irecv];
	*imsg = a->message_ixs;
	*nmsg = a->nix;
}

void history_get_recv_active(const struct history *h, size_t jrecv, const size_t **isend,
			     size_t *nsend)
{
	const struct history_actor *a = &h->receivers[jrecv];
	*isend = a->ind;
	*nsend = a->nz;
}

const size_t *history_recv_counts(const struct history *h, size_t jrecv)
{
	assert(h);
	assert(jrecv < history_recv_count(h));
	const struct history_actor *a = &h->receivers[jrecv];
	return a->nmsg;
}

#endif /* HISTORY_H */
