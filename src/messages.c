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
	
	msgs->tcur = -INFINITY;
	msgs->max_to = -1;
	msgs->max_from = -1;
	msgs->max_nto = 0;
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

struct messages *messages_ref(struct messages * msgs)
{
	assert(msgs);
	refcount_get(&msgs->refcount);
	return msgs;
}

void messages_deinit(struct messages *msgs)
{
	assert(msgs);
	refcount_deinit(&msgs->refcount);
	darray_deinit(&msgs->message_reps);
	darray_deinit(&msgs->recipients);
}

void messages_free(struct messages * msgs)
{
	if (!msgs)
		return;
	
	if (refcount_put(&msgs->refcount, NULL)) {
		messages_deinit(msgs);
		free(msgs);
	}
}

ssize_t iproc_messages_size(struct messages * msgs)
{
	assert(msgs);
	return darray_size(&msgs->message_reps);
}

void iproc_messages_advance_to(struct messages * msgs, double t)
{
	assert(msgs);
	assert(t >= msgs->tcur);
	msgs->tcur = t;
}

void iproc_messages_insert(struct messages * msgs, ssize_t from, ssize_t to)
{
	assert(msgs);
	assert(from >= 0);
	assert(to >= 0);
	iproc_messages_insertm(msgs, from, &to, 1);
}

void
iproc_messages_insertm(struct messages * msgs,
		       ssize_t from, ssize_t *to, ssize_t nto)
{
	assert(msgs);
	assert(from >= 0);
	assert(nto >= 0);
	assert(to || nto == 0);

	double time = msgs->tcur;
	struct darray *message_reps = &msgs->message_reps;
	struct darray *recipients = &msgs->recipients;

	ssize_t ito = darray_size(recipients);
	struct message_rep m = { time, from, ito, nto };
	ssize_t i;

	for (i = 0; i < nto; i++) {
		assert(to[i] >= 0);

		darray_push_back(recipients, to + i);
		if (to[i] > msgs->max_to)
			msgs->max_to = to[i];
	}

	if (from > msgs->max_from)
		msgs->max_from = from;
	if (nto > msgs->max_nto)
		msgs->max_nto = nto;

	darray_push_back(message_reps, &m);
}

ssize_t iproc_messages_max_from(struct messages * msgs)
{
	assert(msgs);
	return msgs->max_from;
}

ssize_t iproc_messages_max_to(struct messages * msgs)
{
	assert(msgs);
	return msgs->max_to;
}

ssize_t iproc_messages_max_nto(struct messages * msgs)
{
	assert(msgs);
	return msgs->max_nto;
}
