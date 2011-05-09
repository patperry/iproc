#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "messages.h"

bool messages_init(struct messages *msgs)
{
	assert(msgs);

	if (!darray_init(&msgs->message_reps, sizeof(struct message_rep)))
		goto fail_message_reps;
	if (!darray_init(&msgs->recipients, sizeof(ssize_t)))
		goto fail_recipients;
	if (!refcount_init(&msgs->refcount))
		goto fail_refcount;

	msgs->tlast = -INFINITY;
	msgs->max_to = -1;
	msgs->max_from = -1;
	msgs->max_nto = 0;
	msgs->to_cached = false;
	return true;

fail_refcount:
	darray_deinit(&msgs->recipients);
fail_recipients:
	darray_deinit(&msgs->message_reps);
fail_message_reps:
	return false;
}

struct messages *messages_alloc()
{
	struct messages *msgs = malloc(sizeof(*msgs));

	if (msgs) {
		if (messages_init(msgs))
			return msgs;
		free(msgs);
	}
	return NULL;
}

struct messages *messages_ref(struct messages *msgs)
{
	assert(msgs);
	if (refcount_get(&msgs->refcount))
		return msgs;
	return NULL;
}

void messages_deinit(struct messages *msgs)
{
	assert(msgs);
	refcount_deinit(&msgs->refcount);
	darray_deinit(&msgs->message_reps);
	darray_deinit(&msgs->recipients);
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

ssize_t messages_size(const struct messages *msgs)
{
	assert(msgs);
	return darray_size(&msgs->message_reps);
}

double messages_tlast(const struct messages *msgs)
{
	assert(msgs);
	return msgs->tlast;
}

struct message *messages_at(const struct messages *msgs, ssize_t i)
{
	assert(msgs);
	assert(0 <= i && i < messages_size(msgs));

	struct message_rep *rep = darray_at(&msgs->message_reps, i);

	if (!msgs->to_cached) {
		ssize_t msg_ito = rep->ito;
		ssize_t *msg_to = darray_at(&msgs->recipients, msg_ito);
		rep->message.to = msg_to;
	}

	return &rep->message;
}

bool messages_add(struct messages *msgs, double time,
		  ssize_t from, ssize_t *to, ssize_t nto, intptr_t attr)
{
	assert(msgs);
	assert(time >= messages_tlast(msgs));
	assert(from >= 0);
	assert(nto >= 0);
	assert(to || nto == 0);

	struct darray *message_reps = &msgs->message_reps;
	struct darray *recipients = &msgs->recipients;

	// reserve space for new recipients and message_reps
	if (!darray_reserve(message_reps, darray_size(message_reps) + 1)
	    || !darray_reserve(recipients, darray_size(recipients) + nto))
		return false;

	ssize_t ito = darray_size(recipients);
	struct message_rep m = { {time, from, NULL, nto, attr}, ito };
	ssize_t i;

	for (i = 0; i < nto; i++) {
		assert(to[i] >= 0);

		darray_push_back(recipients, to + i);	// always succeeds
		if (to[i] > msgs->max_to)
			msgs->max_to = to[i];
	}

	if (from > msgs->max_from)
		msgs->max_from = from;
	if (nto > msgs->max_nto)
		msgs->max_nto = nto;

	darray_push_back(message_reps, &m);	// always succeeds
	msgs->to_cached = false;
	msgs->tlast = time;

	return true;
}

ssize_t messages_max_from(const struct messages *msgs)
{
	assert(msgs);
	return msgs->max_from;
}

ssize_t messages_max_to(const struct messages *msgs)
{
	assert(msgs);
	return msgs->max_to;
}

ssize_t messages_max_nto(const struct messages *msgs)
{
	assert(msgs);
	return msgs->max_nto;
}
