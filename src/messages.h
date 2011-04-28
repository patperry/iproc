#ifndef _MESSAGES_H
#define _MESSAGES_H

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

struct messages_iter {
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


ssize_t messages_size(const struct messages * msgs);
void messages_advance_to(struct messages * msgs, double t);

void messages_insert(struct messages * msgs, ssize_t from, ssize_t to);
void messages_insertm(struct messages * msgs,
			    ssize_t from, ssize_t *to, ssize_t nto);

ssize_t messages_max_from(const struct messages * msgs);
ssize_t messages_max_to(const struct messages * msgs);
ssize_t messages_max_nto(const struct messages * msgs);

struct messages_iter *messages_iter_alloc(struct messages * msgs);
struct messages_iter *messages_iter_ref(struct messages_iter * it);
void messages_iter_free(struct messages_iter * it);

double messages_iter_time(struct messages_iter * it);
ssize_t messages_iter_ntie(struct messages_iter * it);
void messages_iter_select(struct messages_iter * it, ssize_t tie);

ssize_t messages_iter_from(struct messages_iter * it);
ssize_t messages_iter_nto(struct messages_iter * it);
ssize_t *messages_iter_to(struct messages_iter * it);

iproc_history *messages_iter_history(struct messages_iter * it);

void messages_iter_reset(struct messages_iter * it);
bool messages_iter_next(struct messages_iter * it);
bool messages_iter_started(struct messages_iter * it);
bool messages_iter_finished(struct messages_iter * it);

#endif /* _MESSAGES_H */
