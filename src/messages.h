#ifndef _IPROC_MESSAGES_H
#define _IPROC_MESSAGES_H

#include "darray.h"
#include "history.h"
#include "refcount.h"

struct message_rep {
	double time;
	ssize_t from;
	ssize_t ito;
	ssize_t nto;
};

struct messages {
	struct darray message_reps;
	struct darray recipients;
	struct refcount refcount;
	double tcur;
	ssize_t max_from;
	ssize_t max_to;
	ssize_t max_nto;
};

struct message_iter {
	struct messages *messages;
	iproc_history *history;
	ssize_t offset;
	ssize_t ntie;
	struct message_rep *message;
	bool finished;
	struct refcount refcount;
};

bool messages_init(struct messages *msgs);
void messages_deinit(struct messages *msgs);

struct messages *messages_alloc();
struct messages *messages_ref(struct messages *msgs);
void messages_free(struct messages * msgs);


ssize_t iproc_messages_size(struct messages * msgs);
void iproc_messages_advance_to(struct messages * msgs, double t);

void iproc_messages_insert(struct messages * msgs, ssize_t from, ssize_t to);
void iproc_messages_insertm(struct messages * msgs,
			    ssize_t from, ssize_t *to, ssize_t nto);

ssize_t iproc_messages_max_from(struct messages * msgs);
ssize_t iproc_messages_max_to(struct messages * msgs);
ssize_t iproc_messages_max_nto(struct messages * msgs);

struct message_iter *iproc_message_iter_new(struct messages * msgs);
struct message_iter *iproc_message_iter_ref(struct message_iter * it);
void iproc_message_iter_unref(struct message_iter * it);

double iproc_message_iter_time(struct message_iter * it);
ssize_t iproc_message_iter_ntie(struct message_iter * it);
void iproc_message_iter_select(struct message_iter * it, ssize_t tie);

ssize_t iproc_message_iter_from(struct message_iter * it);
ssize_t iproc_message_iter_nto(struct message_iter * it);
ssize_t *iproc_message_iter_to(struct message_iter * it);

iproc_history *iproc_message_iter_history(struct message_iter * it);

void iproc_message_iter_reset(struct message_iter * it);
bool iproc_message_iter_next(struct message_iter * it);
bool iproc_message_iter_started(struct message_iter * it);
bool iproc_message_iter_finished(struct message_iter * it);

#endif /* _IPROC_MESSAGES_H */
