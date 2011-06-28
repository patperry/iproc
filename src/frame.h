#ifndef _FRAME_H
#define _FRAME_H

#include <stdio.h>
#include "actors.h"
#include "messages.h"
#include "history.h"
#include "intmap.h"
#include "pqueue.h"
#include "design.h"

enum frame_event_type {
	DYAD_EVENT_INIT = 1 << 0,
	DYAD_EVENT_MOVE = 1 << 1,
	TRIAD_EVENT_INIT = 1 << 2,
	TRIAD_EVENT_MOVE1 = 1 << 3,
	TRIAD_EVENT_MOVE2 = 1 << 4,
	SEND_VAR_EVENT = 1 << 5,
	RECV_VAR_EVENT = 1 << 6
};

struct message_event_meta {
	ssize_t from;
	ssize_t *to;
	ssize_t nto;
	intptr_t attr;
};

struct dyad_event_meta {
	double msg_time;
	struct dyad msg_dyad;
	intptr_t msg_attr;
	ssize_t intvl;
};

struct triad_event_meta {
	double msg_time1, msg_time2;
	struct dyad msg_dyad1, msg_dyad2;
	intptr_t msg_attr1, msg_attr2;
	ssize_t intvl1, intvl2;
};

struct send_var_event_meta {
	ssize_t item;
	ssize_t index;
	double delta;
};

struct recv_var_event_meta {
	struct dyad item;
	ssize_t index;
	double delta;
};

struct frame_event {
	enum frame_event_type type;
	ssize_t id;
	double time;
	union {
		struct message_event_meta message;
		struct dyad_event_meta dyad;
		struct triad_event_meta triad;
		struct send_var_event_meta send_var;
		struct recv_var_event_meta recv_var;
	} meta;
};

struct frame_handlers {
	uint8_t event_mask;
	void (*handle_event) (void *udata, const struct frame_event *e, struct frame *f);
	void (*handle_clear) (void *udata, struct frame *f);
};

struct frame_observer {
	void *udata;
	struct frame_handlers h;
};

struct frame {
	const struct design *design;
	struct history history;
	struct intmap recv_frames;	// (j, dX[t,i) pairs; dX is a 'struct send_frame'
	struct array vars;
	struct array events;	// events that happen immediately after the current time
	struct pqueue future_events;
	ssize_t next_event_id;
	struct array observers;
	struct matrix send_xt;
	struct refcount refcount;
};

/* dX[t,i] */
struct recv_frame {
	struct frame *frame;
	struct intmap jrecv_dxs;
};

/* create/destroy */
void frame_init(struct frame *f, const struct design *design);
void frame_deinit(struct frame *f);

struct frame *frame_alloc(struct design *design);
struct frame *frame_ref(struct frame *f);
void frame_free(struct frame *f);
void frame_clear(struct frame *f);

/* add a message/advance time  */
void frame_add(struct frame *f, const struct message *msg);
void frame_advance(struct frame *f);
void frame_advance_to(struct frame *f, double t);
double frame_next_change(const struct frame *f);

/* add a future or current event */
ssize_t frame_events_add(struct frame *f, struct frame_event *e);

/* current time/history */
double frame_time(const struct frame *f);
const struct history *frame_history(const struct frame *f);

/* current events */
ssize_t frame_events_count(const struct frame *f);
struct frame_event *frame_events_item(const struct frame *f, ssize_t i);

/* current covariates */
struct vector frame_send_x(struct frame *f, ssize_t isend);
struct vector *frame_recv_dx(const struct frame *f, ssize_t isend,
			     ssize_t jrecv);

/* observers */
void frame_add_observer(struct frame *f, void *udata, const struct frame_handlers *h);
void frame_remove_observer(struct frame *f, void *udata);
			

void frame_send_mul(double alpha, enum trans_op trans,
		    const struct frame *f,
		    const struct vector *x, double beta, struct vector *y);
void frame_send_muls(double alpha, enum trans_op trans,
		     const struct frame *f,
		     const struct svector *x, double beta, struct vector *y);

void frame_recv_mul(double alpha, enum trans_op trans,
		    const struct frame *f, ssize_t isend,
		    const struct vector *x, double beta, struct vector *y);
void frame_recv_muls(double alpha, enum trans_op trans,
		     const struct frame *f, ssize_t isend,
		     const struct svector *x, double beta, struct vector *y);

void frame_recv_dmul(double alpha, enum trans_op trans,
		     const struct frame *f, ssize_t isend,
		     const struct vector *x, double beta, struct svector *y);
void frame_recv_dmuls(double alpha, enum trans_op trans,
		      const struct frame *f, ssize_t isend,
		      const struct svector *x, double beta, struct vector *y);

static inline const struct design *frame_design(const struct frame *f)
{
	assert(f);
	return f->design;
}

/* debugging */
void fprintf_event(FILE * restrict stream, const struct frame_event *e);

#endif /* _FRAME_H */
