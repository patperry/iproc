#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "messages.h"

static void iproc_message_iter_free(struct message_iter * it)
{
	if (it) {
		iproc_history_unref(it->history);
		messages_free(it->messages);
		free(it);
	}
}

struct message_iter *iproc_message_iter_new(struct messages * msgs)
{
	struct message_iter *it = malloc(sizeof(*it));
	if (!it)
		return NULL;

	it->messages = messages_ref(msgs);
	it->history = iproc_history_new();
	refcount_init(&it->refcount);
	iproc_message_iter_reset(it);

	return it;
}

struct message_iter *iproc_message_iter_ref(struct message_iter * it)
{
	if (it) {
		refcount_get(&it->refcount);
	}

	return it;
}

static void iproc_message_iter_release(struct refcount *refcount)
{
	struct message_iter *it =
	    container_of(refcount, struct message_iter, refcount);
	iproc_message_iter_free(it);
}

void iproc_message_iter_unref(struct message_iter * it)
{
	if (it) {
		refcount_put(&it->refcount, iproc_message_iter_release);
	}
}

double iproc_message_iter_time(struct message_iter * it)
{
	assert(iproc_message_iter_started(it));
	assert(!iproc_message_iter_finished(it));
	return it->message->time;
}

ssize_t iproc_message_iter_ntie(struct message_iter * it)
{
	assert(iproc_message_iter_started(it));
	assert(!iproc_message_iter_finished(it));
	return it->ntie;
}

void iproc_message_iter_select(struct message_iter * it, ssize_t tie)
{
	assert(iproc_message_iter_started(it));
	assert(!iproc_message_iter_finished(it));
	assert(tie >= 0);
	assert(tie < iproc_message_iter_ntie(it));

	ssize_t i = it->offset + tie;
	it->message = darray_at(&it->messages->message_reps, i);
}

ssize_t iproc_message_iter_from(struct message_iter * it)
{
	assert(iproc_message_iter_started(it));
	assert(!iproc_message_iter_finished(it));

	return it->message->from;
}

ssize_t iproc_message_iter_nto(struct message_iter * it)
{
	assert(iproc_message_iter_started(it));
	assert(!iproc_message_iter_finished(it));

	return it->message->nto;
}

ssize_t *iproc_message_iter_to(struct message_iter * it)
{
	assert(iproc_message_iter_started(it));
	assert(!iproc_message_iter_finished(it));

	ssize_t ito = it->message->ito;
	ssize_t *to = darray_at(&it->messages->recipients,
				ito);
	return to;
}

void iproc_message_iter_reset(struct message_iter * it)
{
	if (!it)
		return;

	iproc_history_clear(it->history);
	it->message = NULL;
	it->offset = 0;
	it->ntie = 0;
	it->finished = false;
}

bool iproc_message_iter_next(struct message_iter * it)
{
	if (iproc_message_iter_finished(it))
		return false;

	ssize_t offset = it->offset + it->ntie;

	struct darray *message_reps = &it->messages->message_reps;
	struct darray *recipients = &it->messages->recipients;
	iproc_history *history = it->history;
	ssize_t n = darray_size(message_reps);
	bool has_next = offset < n;

	if (has_next) {
		struct message_rep *message = darray_at(message_reps, offset);
		double time = message[0].time;
		ssize_t ntie_max = n - offset;
		ssize_t ntie = 0;

		iproc_history_advance_to(history, time);

		while (ntie < ntie_max && message[ntie].time == time) {
			ssize_t msg_from = message[ntie].from;
			ssize_t msg_ito = message[ntie].ito;
			ssize_t *msg_to = darray_at(recipients, msg_ito);
			ssize_t msg_nto = message[ntie].nto;

			iproc_history_insertm(history, msg_from, msg_to,
					      msg_nto);
			ntie++;
		}
		it->ntie = ntie;
		it->message = message;
	} else {
		it->finished = true;
	}
	it->offset = offset;

	return has_next;
}

bool iproc_message_iter_started(struct message_iter * it)
{
	if (!it)
		return false;

	return it->message != NULL;
}

bool iproc_message_iter_finished(struct message_iter * it)
{
	if (!it)
		return true;

	return it->finished;
}

iproc_history *iproc_message_iter_history(struct message_iter * it)
{
	if (!it)
		return NULL;

	return it->history;
}
