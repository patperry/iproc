#ifndef _FRAME_H
#define _FRAME_H

#include "actors.h"
#include "messages.h"
#include "intmap.h"
#include "pqueue.h"
#include "design.h"

struct frame {
	struct pqueue dyad_var_diffs;
	struct intmap send_frames; // (j, dX[t,i) pairs; dX is a 'struct send_frame'
	struct design *design;
	double time;
};

/* dX[t,i] */
struct send_frame {
	struct frame *frame;
	struct intmap jrecv_dxs; // (j, dx[t,i,j]) pairs; dx is a 'struct svector'
};

struct dyad_var_diff {
	double time;	
	double delta;
	ssize_t index;	
	ssize_t isend;
	ssize_t jrecv;
};

struct frame_diff {
	struct darray dyad_var_diffs;
};

/* create/destroy */
bool frame_init(struct frame *f, struct design *design);
void frame_deinit(struct frame *f);

bool frame_add_dyad_event(struct frame *f, const struct dyad_var_diff *delta);
bool frame_reserve_dyad_events(struct frame *f, ssize_t nadd);

/* record a message event */
bool frame_insert(struct frame *f, const struct message *msg);

/* advance time */
bool frame_advance_to(struct frame *f, double t, struct frame_diff *diff);

void frame_clear(struct frame *f);

/* time of the next change */
double frame_time(const struct frame *f);
double frame_next_update(const struct frame *f);



bool frame_mul(double alpha, enum trans_op trans,
	       const struct frame *f, ssize_t isend,
	       const struct vector *x, double beta, struct vector *y);
bool frame_muls(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct svector *x, double beta, struct vector *y);

bool frame_dmul(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct vector *x, double beta, struct svector *y);
bool frame_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f, ssize_t isend,
		 const struct svector *x, double beta, struct svector *y);



#endif /* _FRAME_H */