#include "port.h"
#include <assert.h>
#include <math.h>
#include "ieee754.h"
#include "dyad-queue.h"

static int dyad_queue_event_rcompare(const void *x, const void *y)
{
	const struct dyad_queue_event *e = x;
	const struct dyad_queue_event *f = y;
	return double_rcompare(&e->tnext, &f->tnext);
}

bool dyad_queue_init(struct dyad_queue *queue, const struct vector *intervals)
{
	assert(queue);
	assert(intervals);
	assert(vector_empty(intervals) || vector_get_front(intervals) > 0.0);
	assert(vector_empty(intervals) || vector_get_back(intervals) < INFINITY);
#ifndef NDEBUG
	ssize_t i;
	for (i = 1; i < vector_dim(intervals); i++) {
		assert(vector_get(intervals, i - 1) < vector_get(intervals, i));
	}
#endif
	
	if (!pqueue_init(&queue->events, dyad_queue_event_rcompare,
			 sizeof(struct dyad_queue_event)))
		goto fail_events;
	
	queue->intervals = intervals;
	queue->next_id = 0;
	return true;
	
fail_events:
	return false;
}

void dyad_queue_deinit(struct dyad_queue *queue)
{
	assert(queue);
	pqueue_deinit(&queue->events);
}

void dyad_queue_clear(struct dyad_queue *queue)
{
	assert(queue);
	pqueue_clear(&queue->events);
	queue->next_id = 0;
}

bool dyad_queue_empty(const struct dyad_queue *queue)
{
	assert(queue);
	return pqueue_empty(&queue->events);
}

ssize_t dyad_queue_size(const struct dyad_queue *queue)
{
	assert(queue);
	return pqueue_size(&queue->events);
}

double dyad_queue_next_update(const struct dyad_queue *queue)
{
	assert(queue);
	
	if (dyad_queue_empty(queue))
		return INFINITY;
	
	const struct dyad_queue_event *top = pqueue_top(&queue->events);
	return top->tnext;
}


const struct dyad_event *dyad_queue_top(const struct dyad_queue *queue)
{
	assert(queue);
	assert(!dyad_queue_empty(queue));
	
	const struct dyad_queue_event *top = pqueue_top(&queue->events);
	return &top->event;
}

bool dyad_queue_push(struct dyad_queue *queue, const struct message *msg)
{
	assert(queue);
	assert(msg);
	assert(isfinite(msg->time));	
	assert(!isnan(msg->time));
	assert(msg->nto >= 0);
	assert(queue->next_id <= SSIZE_MAX - msg->nto);
	
	ssize_t ito, nto = msg->nto;
	struct dyad_queue_event e;
	
	e.tnext = msg->time;
	e.event.type = DYAD_EVENT_INIT;
	e.event.time = msg->time;
	e.event.isend = msg->from;
	e.event.attr = msg->attr;
	e.event.intvl = 0;
	
	for (ito = 0; ito < nto; ito++) {
		e.event.id = queue->next_id++;		
		e.event.jrecv = msg->to[ito];
		
		if (!pqueue_push(&queue->events, &e))
			return false;
	}
	
	return true;
}

void dyad_queue_pop(struct dyad_queue *queue)
{
	assert(queue);
	
	struct dyad_queue_event *e = pqueue_top(&queue->events);
	double t0, dt;
	
	if (e->event.intvl == vector_dim(queue->intervals)) {
		pqueue_pop(&queue->events);
	} else {
		t0 = e->event.time;
		dt = vector_get(queue->intervals, e->event.intvl);
		
		e->tnext = t0 + dt;
		e->event.type = DYAD_EVENT_MOVE;
		e->event.intvl++;
		
		pqueue_update_top(&queue->events);
	}
}
