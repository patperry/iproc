#ifndef _MESSAGES_H
#define _MESSAGES_H

#include <stddef.h>
#include <stdint.h>

struct message {
	double time;
	size_t from;
	size_t *to;
	size_t nto;
	intptr_t attr;
};

struct message_rep {
	struct message message;
	size_t ito;
};

struct messages {
	struct message_rep *reps;
	size_t nsend, nsend_max;
	size_t *recv;
	size_t nrecv, nrecv_max;
	double tlast;
	size_t max_from;
	size_t max_to;
	size_t max_nto;
	int to_cached;
};

struct messages_iter {
	struct messages *messages;
	ptrdiff_t offset;
	size_t ntie;
	struct message_rep *message_rep;
};

#define MESSAGES_TIME(it) ((it).message_rep->message.time)
#define MESSAGES_COUNT(it) ((it).ntie)
#define MESSAGES_VAL(it, i) (&(it).message_rep[i].message)
#define MESSAGES_FOREACH(it, msgs) \
	for ((it) = messages_iter_make(msgs); messages_iter_advance(&(it));)

void messages_init(struct messages *msgs);
void messages_deinit(struct messages *msgs);

size_t messages_count(const struct messages *msgs);
size_t messages_recv_count(const struct messages *msgs);
double messages_tlast(const struct messages *msgs);
void messages_add(struct messages *msgs, double time,
		  size_t from, size_t *to, size_t nto, intptr_t attr);

struct message *messages_at(const struct messages *msgs, size_t i);

size_t messages_max_from(const struct messages *msgs);
size_t messages_max_to(const struct messages *msgs);
size_t messages_max_nto(const struct messages *msgs);

struct messages_iter messages_iter_make(const struct messages *msgs);
void messages_iter_reset(struct messages_iter *it);
int messages_iter_advance(struct messages_iter *it);

#endif /* _MESSAGES_H */