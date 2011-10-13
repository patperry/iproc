#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "xalloc.h"
#include "messages.h"

void messages_init(struct messages *msgs)
{
	assert(msgs);

	array_init(&msgs->message_reps, sizeof(struct message_rep));
	array_init(&msgs->recipients, sizeof(size_t));
	refcount_init(&msgs->refcount);

	msgs->nrecv = 0;
	msgs->tlast = -INFINITY;
	msgs->max_to = 0;
	msgs->max_from = 0;
	msgs->max_nto = 0;
	msgs->to_cached = false;
}

struct messages *messages_alloc()
{
	struct messages *msgs = xcalloc(1, sizeof(*msgs));
	messages_init(msgs);
	return msgs;
}

struct messages *messages_ref(struct messages *msgs)
{
	assert(msgs);
	refcount_get(&msgs->refcount);
	return msgs;
}

void messages_deinit(struct messages *msgs)
{
	assert(msgs);
	refcount_deinit(&msgs->refcount);
	array_deinit(&msgs->message_reps);
	array_deinit(&msgs->recipients);
}

void messages_free(struct messages *msgs)
{
	if (!msgs)
		return;

	if (refcount_put(&msgs->refcount, NULL)) {
		refcount_get(&msgs->refcount);
		messages_deinit(msgs);
		free(msgs);
	}
}

size_t messages_count(const struct messages *msgs)
{
	assert(msgs);
	return array_count(&msgs->message_reps);
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

	struct message_rep *rep = array_item(&msgs->message_reps, i);

	if (!msgs->to_cached) {
		size_t msg_ito = rep->ito;
		size_t *msg_to = array_item(&msgs->recipients, msg_ito);
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

	struct array *message_reps = &msgs->message_reps;
	struct array *recipients = &msgs->recipients;

	size_t ito = array_count(recipients);
	struct message_rep m = { {time, from, NULL, nto, attr}, ito };
	size_t i;

	for (i = 0; i < nto; i++) {
		array_add(recipients, to + i);	// always succeeds
		if (to[i] > msgs->max_to)
			msgs->max_to = to[i];
	}

	if (from > msgs->max_from)
		msgs->max_from = from;
	if (nto > msgs->max_nto)
		msgs->max_nto = nto;

	array_add(message_reps, &m);	// always succeeds
	msgs->to_cached = false;
	msgs->tlast = time;
	msgs->nrecv += nto;
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

bool messages_iter_advance(struct messages_iter *it)
{
	assert(it);
	assert(it->offset < (ptrdiff_t)messages_count(it->messages));

	size_t offset = it->offset + it->ntie;

	struct array *message_reps = &it->messages->message_reps;
	struct array *recipients = &it->messages->recipients;
	size_t n = array_count(message_reps);
	bool has_next = offset < n;

	if (has_next) {
		struct message_rep *message_rep =
		    array_item(message_reps, offset);
		double time = message_rep[0].message.time;
		size_t ntie_max = n - offset;
		size_t ntie = 0;

		do {
			if (!it->messages->to_cached) {
				size_t msg_ito = message_rep[ntie].ito;
				size_t *msg_to =
				    array_item(recipients, msg_ito);
				message_rep[ntie].message.to = msg_to;
			}

			ntie++;
		} while (ntie < ntie_max
			 && message_rep[ntie].message.time == time);

		it->ntie = ntie;
		it->message_rep = message_rep;
	} else {
		it->messages->to_cached = true;
	}
	it->offset = offset;

	return has_next;
}
