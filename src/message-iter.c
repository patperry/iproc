#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "messages.h"



struct messages_iter messages_iter(struct messages * msgs)
{
	struct messages_iter it;

	it.messages = messages_ref(msgs);
	it.history = iproc_history_new();
	messages_iter_reset(&it);

	return it;
}

void messages_iter_deinit(struct messages_iter *it)
{
	iproc_history_unref(it->history);
	messages_free(it->messages);
}

double messages_iter_current_time(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));
	return it->message_rep->message.time;
}

ssize_t messages_iter_ntie(struct messages_iter * it)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));
	return it->ntie;
}

struct message *messages_iter_current(struct messages_iter *it, ssize_t itie)
{
	assert(messages_iter_started(it));
	assert(!messages_iter_finished(it));
	assert(itie >= 0);
	assert(itie < messages_iter_ntie(it));

	ssize_t i = it->offset + itie;
	it->message_rep = darray_at(&it->messages->message_reps, i);
	return &it->message_rep->message;

}

void messages_iter_reset(struct messages_iter * it)
{
	if (!it)
		return;

	iproc_history_clear(it->history);
	it->message_rep = NULL;
	it->offset = 0;
	it->ntie = 0;
	it->finished = false;
}

bool messages_iter_advance(struct messages_iter * it)
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
		struct message_rep *message_rep = darray_at(message_reps, offset);
		double time = message_rep[0].message.time;
		ssize_t ntie_max = n - offset;
		ssize_t ntie = 0;

		iproc_history_advance_to(history, time);

		do {
			if (!it->messages->to_cached) {
				ssize_t msg_ito = message_rep[ntie].ito;
				ssize_t *msg_to = darray_at(recipients, msg_ito);
				message_rep[ntie].message.to = msg_to;
			}
			
			/* deprecated */
			ssize_t msg_from = message_rep[ntie].message.from;
			ssize_t *msg_to = message_rep[ntie].message.to;
			ssize_t msg_nto = message_rep[ntie].message.nto;
			intptr_t attr = message_rep[ntie].message.attr;
			iproc_history_insertm(history, msg_from, msg_to,
					      msg_nto, attr);
			
			/* not deprecated */
			ntie++;
		} while (ntie < ntie_max && message_rep[ntie].message.time == time);
		
		it->ntie = ntie;
		it->message_rep = message_rep;
	} else {
		it->finished = true;
		it->messages->to_cached = true;
	}
	it->offset = offset;

	return has_next;
}

bool messages_iter_started(struct messages_iter * it)
{
	if (!it)
		return false;

	return it->message_rep != NULL;
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
