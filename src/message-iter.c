#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "messages.h"

struct messages_iter messages_iter(struct messages *msgs)
{
	assert(msgs);

	struct messages_iter it;

	it.messages = msgs;
	messages_iter_reset(&it);

	return it;
}

double messages_iter_current_time(struct messages_iter *it)
{
	assert(it);
	assert(0 <= it->offset && it->offset < messages_size(it->messages));

	return it->message_rep->message.time;
}

ssize_t messages_iter_ntie(struct messages_iter *it)
{
	assert(it);
	assert(0 <= it->offset && it->offset < messages_size(it->messages));

	return it->ntie;
}

struct message *messages_iter_current(struct messages_iter *it, ssize_t itie)
{
	assert(it);
	assert(0 <= it->offset && it->offset < messages_size(it->messages));
	assert(0 <= itie && itie < messages_iter_ntie(it));

	ssize_t i = it->offset + itie;
	it->message_rep = darray_at(&it->messages->message_reps, i);
	return &it->message_rep->message;

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
	assert(it->offset < messages_size(it->messages));

	ssize_t offset = it->offset + it->ntie;

	struct darray *message_reps = &it->messages->message_reps;
	struct darray *recipients = &it->messages->recipients;
	ssize_t n = darray_size(message_reps);
	bool has_next = offset < n;

	if (has_next) {
		struct message_rep *message_rep =
		    darray_at(message_reps, offset);
		double time = message_rep[0].message.time;
		ssize_t ntie_max = n - offset;
		ssize_t ntie = 0;

		do {
			if (!it->messages->to_cached) {
				ssize_t msg_ito = message_rep[ntie].ito;
				ssize_t *msg_to =
				    darray_at(recipients, msg_ito);
				message_rep[ntie].message.to = msg_to;
			}

			/* not deprecated */
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
