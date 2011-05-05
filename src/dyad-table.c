#include "port.h"
#include <assert.h>
#include <math.h>
#include "ieee754.h"
#include "dyad-table.h"

static int dyad_queue_event_rcompare(const void *x, const void *y)
{
	const struct dyad_queue_event *e = x;
	const struct dyad_queue_event *f = y;
	return double_rcompare(&e->tnext, &f->tnext);
}

static bool dyad_queue_callback_equals(const void *x, const void *y)
{
	const struct dyad_queue_callback *c = x;
	const struct dyad_queue_callback *d = y;
	
	return (c->udata == d->udata
		&& c->callback == d->callback
		&& c->type == d->type);
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
	
	if (!darray_init(&queue->callbacks, sizeof(struct dyad_queue_callback)))
		goto fail_callbacks;
	
	queue->intervals = intervals;
	queue->time = -INFINITY;
	queue->next_id = 0;
	return true;
	
fail_callbacks:
	pqueue_deinit(&queue->events);
fail_events:
	return false;
}

void dyad_queue_deinit(struct dyad_queue *queue)
{
	assert(queue);
	darray_deinit(&queue->callbacks);
	pqueue_deinit(&queue->events);
}

void dyad_queue_clear(struct dyad_queue *queue)
{
	assert(queue);
	pqueue_clear(&queue->events);
	queue->time = -INFINITY;
	queue->next_id = 0;
}

bool dyad_queue_add_callback(struct dyad_queue *queue,
			     enum dyad_event_type type,
			     dyad_event_fn callback,
			     void *udata)
{
	assert(queue);
	assert(callback);
	
	struct dyad_queue_callback *c;
	
	if ((c = darray_push_back(&queue->callbacks, NULL))) {
		c->type = type;
		c->callback = callback;
		c->udata = udata;
		return true;
	}
	return false;
}

void dyad_queue_remove_callback(struct dyad_queue *queue,
				enum dyad_event_type type,
				dyad_event_fn callback,
				void *udata)
{
	assert(queue);
	assert(callback);
	
	struct dyad_queue_callback key = { type, callback, udata };
	ssize_t index = darray_find_last_index(&queue->callbacks, &key,
					       dyad_queue_callback_equals);
	
	if (index >= 0) {
		darray_erase(&queue->callbacks, index);
	}
}

void dyad_queue_clear_callbacks(struct dyad_queue *queue)
{
	assert(queue);
	darray_clear(&queue->callbacks);
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

double dyad_queue_time(const struct dyad_queue *queue)
{
	assert(queue);
	return queue->time;
}

double dyad_queue_next(const struct dyad_queue *queue)
{
	assert(queue);
	
	if (pqueue_empty(&queue->events))
		return INFINITY;
	
	const struct dyad_queue_event *top = pqueue_top(&queue->events);
	return top->tnext;
}

bool dyad_queue_insert(struct dyad_queue *queue, const struct message *msg)
{
	assert(queue);
	assert(msg);
	assert(isfinite(msg->time));	
	assert(msg->time == dyad_queue_time(queue));
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

bool dyad_queue_advance_to(struct dyad_queue *queue, double time)
{
	assert(queue);
	assert(time >= dyad_queue_time(queue));
	
	ssize_t i, n;
	struct dyad_queue_event *e;
	const struct dyad_queue_callback *c;
	double t0, dt;
	bool ok;
	
	if (time == dyad_queue_time(queue))
		return true;
	
	n = darray_size(&queue->callbacks);
	
	while (dyad_queue_next(queue) < time) {
		e = pqueue_top(&queue->events);
		
		// notify all listeners
		for (i = 0; i < n; i++) {
			c = darray_at(&queue->callbacks, i);
			if (c->type == e->event.type) {
				ok = c->callback(&e->event, c->udata);
				if (!ok)
					return false;
			}
		}
		
		// update event
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
	queue->time = time;
	
	return true;
}
