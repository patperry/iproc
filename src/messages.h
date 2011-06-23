#ifndef _MESSAGES_H
#define _MESSAGES_H

#include "array.h"
#include "history.h"
#include "refcount.h"

struct message {
	double time;
	ssize_t from;
	ssize_t *to;
	ssize_t nto;
	intptr_t attr;
};

struct message_rep {
	struct message message;
	ssize_t ito;
};

struct messages {
	struct array message_reps;
	struct array recipients;
	struct refcount refcount;
	ssize_t nrecv;
	double tlast;
	ssize_t max_from;
	ssize_t max_to;
	ssize_t max_nto;
	bool to_cached;
};

struct messages_iter {
	struct messages *messages;
	ssize_t offset;
	ssize_t ntie;
	struct message_rep *message_rep;
};

#define MESSAGES_TIME(it) ((it).message_rep->message.time)
#define MESSAGES_COUNT(it) ((it).ntie)
#define MESSAGES_VAL(it, i) (&(it).message_rep[i].message)
#define MESSAGES_FOREACH(it, msgs) \
	for ((it) = messages_iter_make(msgs); messages_iter_advance(&(it));)

void messages_init(struct messages *msgs);
void messages_deinit(struct messages *msgs);

struct messages *messages_alloc();
struct messages *messages_ref(struct messages *msgs);
void messages_free(struct messages *msgs);

ssize_t messages_count(const struct messages *msgs);
ssize_t messages_recv_count(const struct messages *msgs);
double messages_tlast(const struct messages *msgs);
bool messages_add(struct messages *msgs, double time,
		  ssize_t from, ssize_t *to, ssize_t nto, intptr_t attr);

struct message *messages_at(const struct messages *msgs, ssize_t i);

ssize_t messages_max_from(const struct messages *msgs);
ssize_t messages_max_to(const struct messages *msgs);
ssize_t messages_max_nto(const struct messages *msgs);

struct messages_iter messages_iter_make(const struct messages *msgs);
void messages_iter_reset(struct messages_iter *it);
bool messages_iter_advance(struct messages_iter *it);

#endif /* _MESSAGES_H */
