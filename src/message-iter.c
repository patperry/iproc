#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "messages.h"


struct messages_iter *messages_iter_alloc(struct messages * msgs)
{
	struct messages_iter *it = malloc(sizeof(*it));
	if (!it)
		return NULL;

	it->messages = messages_ref(msgs);
	it->history = iproc_history_new();
	messages_iter_reset(it);

	return it;
}

void message_iter_deinit(struct messages_iter *it)
{
	iproc_history_unref(it->history);
	messages_free(it->messages);
}

void messages_iter_free(struct messages_iter * it)
{
	if (it) {
		message_iter_deinit(it);
		free(it);
	}
}

double messages_iter_time(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));
	return it->message->time;
}

ssize_t messages_iter_ntie(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));
	return it->ntie;
}

void messages_iter_select(struct messages_iter * it, ssize_t tie)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));
	assert(tie >= 0);
	assert(tie < messages_iter_ntie(it));

	ssize_t i = it->offset + tie;
	it->message = darray_at(&it->messages->message_reps, i);
}

ssize_t messages_iter_from(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));

	return it->message->from;
}

ssize_t messages_iter_nto(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));

	return it->message->nto;
}

ssize_t *messages_iter_to(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));

	ssize_t ito = it->message->ito;
	ssize_t *to = darray_at(&it->messages->recipients,
				ito);
	return to;
}

void messages_iter_reset(struct messages_iter * it)
{
	if (!it)
		return;

	iproc_history_clear(it->history);
	it->message = NULL;
	it->offset = 0;
	it->ntie = 0;
	it->finished = false;
}

bool messages_iter_next(struct messages_iter * it)
{
	if (messages_iter_finished(it))
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

bool messages_iter_started(struct messages_iter * it)
{
	if (!it)
		return false;

	return it->message != NULL;
}

bool messages_iter_finished(struct messages_iter * it)
{
	if (!it)
		return true;

	return it->finished;
}

iproc_history *messages_iter_history(struct messages_iter * it)
{
	if (!it)
		return NULL;

	return it->history;
}
