#ifndef _DYAD_QUEUE_H
#define _DYAD_QUEUE_H

#include "history.h"
#include "messages.h"
#include "pqueue.h"
#include "vector.h"

enum dyad_event_type {
	DYAD_EVENT_INIT = 1 << 0,
	DYAD_EVENT_MOVE = 1 << 1
};

struct dyad {
	ssize_t isend;
	ssize_t jrecv;
};

struct dyad_event {
	ssize_t id;
	enum dyad_event_type type;
	struct dyad dyad;
	double time;
	intptr_t attr;
	ssize_t intvl;
};

struct dyad_queue_event {
	double tnext;
	struct dyad_event event;
};

struct dyad_queue {
	struct pqueue events;
	const struct vector *intervals;
	ssize_t next_id;
};

enum triad_event_type {
	TRIAD_EVENT_INIT = 1 << 0,
	TRIAD_EVENT_MOVE1 = 1 << 1,
	TRIAD_EVENT_MOVE2 = 1 << 2,
};

struct triad_event {
	ssize_t id;
	enum triad_event_type type;
	double time1, time2;
	ssize_t isend1, isend2;
	ssize_t jrecv1, jrecv2;
	intptr_t attr1, attr2;
	ssize_t intvl1, intvl2;
};

void dyad_queue_init(struct dyad_queue *queue, const struct vector *intervals);
void dyad_queue_deinit(struct dyad_queue *queue);
void dyad_queue_clear(struct dyad_queue *queue);

bool dyad_queue_empty(const struct dyad_queue *queue);
ssize_t dyad_queue_size(const struct dyad_queue *queue);

double dyad_queue_next_update(const struct dyad_queue *queue);
const struct dyad_event *dyad_queue_top(const struct dyad_queue *queue);

void dyad_queue_push(struct dyad_queue *queue, const struct message *msg);
void dyad_queue_pop(struct dyad_queue *queue);

#endif /* _DYAD_QUEUE_H */
