#ifndef _IPROC_MESSAGES_H
#define _IPROC_MESSAGES_H

#include <stdbool.h>
#include <stdint.h>
#include "darray.h"
#include "history.h"
#include "refcount.h"

typedef struct _iproc_message iproc_message;
typedef struct _iproc_messages iproc_messages;
typedef struct _iproc_message_iter iproc_message_iter;

struct _iproc_message {
	double time;
	ssize_t from;
	ssize_t ito;
	ssize_t nto;
};

struct _iproc_messages {
	double tcur;
	struct darray array;
	struct darray recipients;
	struct refcount refcount;
	ssize_t max_from;
	ssize_t max_to;
	ssize_t max_nto;
};

struct _iproc_message_iter {
	iproc_messages *messages;
	iproc_history *history;
	ssize_t offset;
	ssize_t ntie;
	iproc_message *message;
	bool finished;
	struct refcount refcount;
};

iproc_messages *iproc_messages_new();
iproc_messages *iproc_messages_ref(iproc_messages * msgs);
void iproc_messages_unref(iproc_messages * msgs);

ssize_t iproc_messages_size(iproc_messages * msgs);
void iproc_messages_advance_to(iproc_messages * msgs, double t);

void iproc_messages_insert(iproc_messages * msgs, ssize_t from, ssize_t to);
void iproc_messages_insertm(iproc_messages * msgs,
			    ssize_t from, ssize_t *to, ssize_t nto);

ssize_t iproc_messages_max_from(iproc_messages * msgs);
ssize_t iproc_messages_max_to(iproc_messages * msgs);
ssize_t iproc_messages_max_nto(iproc_messages * msgs);

iproc_message_iter *iproc_message_iter_new(iproc_messages * msgs);
iproc_message_iter *iproc_message_iter_ref(iproc_message_iter * it);
void iproc_message_iter_unref(iproc_message_iter * it);

double iproc_message_iter_time(iproc_message_iter * it);
ssize_t iproc_message_iter_ntie(iproc_message_iter * it);
void iproc_message_iter_select(iproc_message_iter * it, ssize_t tie);

ssize_t iproc_message_iter_from(iproc_message_iter * it);
ssize_t iproc_message_iter_nto(iproc_message_iter * it);
ssize_t *iproc_message_iter_to(iproc_message_iter * it);

iproc_history *iproc_message_iter_history(iproc_message_iter * it);

void iproc_message_iter_reset(iproc_message_iter * it);
bool iproc_message_iter_next(iproc_message_iter * it);
bool iproc_message_iter_started(iproc_message_iter * it);
bool iproc_message_iter_finished(iproc_message_iter * it);

#endif /* _IPROC_MESSAGES_H */
