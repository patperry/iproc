#ifndef _FRAME_H
#define _FRAME_H

#include "actors.h"
#include "messages.h"
#include "history.h"
#include "intmap.h"
#include "pqueue.h"
#include "design.h"

struct frame {
	struct history history;
	struct dyad_queue dyad_queue;
	struct intmap send_frames; // (j, dX[t,i) pairs; dX is a 'struct send_frame'
	struct design *design;
	double time;
	ssize_t isend;
	struct intmap jrecv_dxs; // dx is a 'struct svector *'
	struct darray dxs;
	ssize_t next_dx;
};

/* dX[t,i] */
struct send_frame {
	struct frame *frame;
	struct intmap jrecv_dxs;
};

/* create/destroy */
bool frame_init(struct frame *f, struct design *design);
void frame_deinit(struct frame *f);

/* record a message event */
bool frame_insert(struct frame *f, const struct message *msg);

/* advance time */
bool frame_advance_to(struct frame *f, double t);

ssize_t frame_sender(const struct frame *f);
bool frame_set_sender(struct frame *f, ssize_t isend);


void frame_clear(struct frame *f);

/* time of the next change */
double frame_time(const struct frame *f);
double frame_next_update(const struct frame *f);
const struct history *frame_history(const struct frame *f);


struct svector *frame_dx(struct frame *f, ssize_t isend, ssize_t jrecv);
struct svector *DEPRECATED_frame_dx(struct frame *f, ssize_t jrecv);


bool frame_mul(double alpha, enum trans_op trans,
	       const struct frame *f,
	       const struct vector *x, double beta, struct vector *y);
bool frame_muls(double alpha, enum trans_op trans,
		const struct frame *f,
		const struct svector *x, double beta, struct vector *y);

bool frame_dmul(double alpha, enum trans_op trans,
		const struct frame *f,
		const struct vector *x, double beta, struct svector *y);
bool frame_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f,
		 const struct svector *x, double beta, struct svector *y);



#endif /* _FRAME_H */