#ifndef HISTORY_H
#define HISTORY_H

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include "uintset.h"
#include "version.h"


struct message {
	double time;
	size_t from;
	size_t *to;
	size_t nto;
	intptr_t attr;
};

struct history_actor {
	struct uintset msgs;
	struct uintset alters;
	double *wts;
};


struct history {
	size_t nsend, nrecv;
	struct message *msgs;
	size_t nmsg, nmsg_max;
	size_t *to;
	size_t nto, nto_max;
	struct history_actor *sends;
	struct history_actor *recvs;
	double tcur;
	size_t ncur;
	struct version version;
};


/* create/destroy/clear */
void history_init(struct history *h, size_t nsend, size_t nrecv);
void history_deinit(struct history *h);
void history_clear(struct history *h); // delete all messages, reset time

/* properties */
static inline size_t history_send_count(const struct history *h) { return h->nsend; }
static inline size_t history_recv_count(const struct history *h) { return h->nrecv; }

/* access past messages */
static inline size_t history_count(const struct history *h) { return h->ncur; }


static inline struct message *history_item(const struct history *h, size_t i)
{
	assert(i < history_count(h));
	assert(h->msgs[i].time < h->tcur);
	return &h->msgs[i];
}

/* modify time */
static inline double history_time(const struct history *h) { return h->tcur; }
void history_reset(struct history *h); // reset time to -INFINITY
void history_advance(struct history *h, double time);	// advance time
void history_set_time(struct history *h, double time);


/* add messages */
void history_add(struct history *h, size_t from, const size_t *to, size_t nto,
		 intptr_t attr);

/* message array access */
static inline void history_get_messages(const struct history *h, const struct message **msgs, size_t *n)
{
	*msgs = h->msgs;
	*n = h->nmsg;
}

/* version (inremented on history_reset or history_clear) */
static inline struct version *history_version(const struct history *h)
{
	return &((struct history *)h)->version;
}


/* actor histories */
static inline struct history_actor *history_send(const struct history *h, size_t i)
{
	assert(i < history_send_count(h));
	return &h->sends[i];
}

static inline struct history_actor *history_recv(const struct history *h, size_t j)
{
	assert(j < history_recv_count(h));
	return &h->recvs[j];
}

static inline void history_actor_get_msgs(const struct history_actor *ha,
					  const size_t **ind, size_t *len)
{
	uintset_get_vals(&ha->msgs, ind, len);
}

static inline void history_actor_get_alters(const struct history_actor *ha,
					    const double **wts, const size_t **ind,
					    size_t *len)
{
	uintset_get_vals(&ha->alters, ind, len);
	*wts = ha->wts;
}


#endif /* HISTORY_H */
