#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include "coreutil.h"
#include "xalloc.h"

#include "messages.h"


static void messages_grow_reps(struct messages *msgs, size_t delta)
{
	size_t nmax = array_grow(msgs->nsend, msgs->nsend_max, delta, SIZE_MAX);
	if (nmax > msgs->nsend_max) {
		msgs->reps = xrealloc(msgs->reps, nmax * sizeof(msgs->reps[0]));
		msgs->nsend_max = nmax;
	}
}

static void messages_grow_recv(struct messages *msgs, size_t delta)
{
	size_t nmax = array_grow(msgs->nrecv, msgs->nrecv_max, delta, SIZE_MAX);
	if (nmax > msgs->nrecv_max) {
		msgs->recv = xrealloc(msgs->recv, nmax * sizeof(msgs->recv[0]));
		msgs->nrecv_max = nmax;
	}
}

void messages_init(struct messages *msgs)
{
	assert(msgs);

	msgs->reps = NULL;
	msgs->nsend = 0;
	msgs->nsend_max = 0;
	msgs->recv = NULL;
	msgs->nrecv = 0;
	msgs->nrecv_max = 0;
	msgs->tlast = -INFINITY;
	msgs->max_to = 0;
	msgs->max_from = 0;
	msgs->max_nto = 0;
	msgs->to_cached = 0;
}

void messages_deinit(struct messages *msgs)
{
	assert(msgs);
	free(msgs->reps);
	free(msgs->recv);
}

size_t messages_count(const struct messages *msgs)
{
	assert(msgs);
	return msgs->nsend;
}

size_t messages_recv_count(const struct messages *msgs)
{
	assert(msgs);
	return msgs->nrecv;
}

double messages_tlast(const struct messages *msgs)
{
	assert(msgs);
	return msgs->tlast;
}

struct message *messages_at(const struct messages *msgs, size_t i)
{
	assert(msgs);
	assert(i < messages_count(msgs));

	struct message_rep *rep = &msgs->reps[i];

	if (!msgs->to_cached) {
		size_t msg_ito = rep->ito;
		size_t *msg_to = &msgs->recv[msg_ito];
		rep->message.to = msg_to;
	}

	return &rep->message;
}

void messages_add(struct messages *msgs, double time,
		  size_t from, size_t *to, size_t nto, intptr_t attr)
{
	assert(msgs);
	assert(time >= messages_tlast(msgs));
	assert(to || nto == 0);

	size_t ito = msgs->nrecv;
	struct message_rep m = { {time, from, NULL, nto, attr}, ito };
	size_t i;

	msgs->tlast = time;

	if (from > msgs->max_from)
		msgs->max_from = from;
	if (nto > msgs->max_nto)
		msgs->max_nto = nto;

	messages_grow_recv(msgs, nto);
	for (i = 0; i < nto; i++, ito++) {
		if (to[i] > msgs->max_to)
			msgs->max_to = to[i];
		msgs->recv[ito] = to[i];
	}
	msgs->nrecv = ito;

	messages_grow_reps(msgs, 1);
	msgs->reps[msgs->nsend++] = m;
	msgs->to_cached = 0;
}

size_t messages_max_from(const struct messages *msgs)
{
	assert(msgs);
	return msgs->max_from;
}

size_t messages_max_to(const struct messages *msgs)
{
	assert(msgs);
	return msgs->max_to;
}

size_t messages_max_nto(const struct messages *msgs)
{
	assert(msgs);
	return msgs->max_nto;
}

struct messages_iter messages_iter_make(const struct messages *msgs)
{
	assert(msgs);

	struct messages_iter it;

	it.messages = (struct messages *)msgs;
	messages_iter_reset(&it);

	return it;
}

void messages_iter_reset(struct messages_iter *it)
{
	assert(it);

	it->message_rep = NULL;
	it->offset = -1;
	it->ntie = 1;
}

int messages_iter_advance(struct messages_iter *it)
{
	assert(it);
	assert(it->offset < (ptrdiff_t)messages_count(it->messages));

	size_t offset = it->offset + it->ntie;

	struct messages *msgs = it->messages;
	size_t n = msgs->nsend;
	int has_next = offset < n;

	if (has_next) {
		struct message_rep *message_rep = &msgs->reps[offset];
		double time = message_rep[0].message.time;
		size_t ntie_max = n - offset;
		size_t ntie = 0;

		do {
			if (!it->messages->to_cached) {
				size_t msg_ito = message_rep[ntie].ito;
				size_t *msg_to = &msgs->recv[msg_ito];
				message_rep[ntie].message.to = msg_to;
			}

			ntie++;
		} while (ntie < ntie_max
			 && message_rep[ntie].message.time == time);

		it->ntie = ntie;
		it->message_rep = message_rep;
	} else {
		it->messages->to_cached = 1;
	}
	it->offset = offset;

	return has_next;
}
