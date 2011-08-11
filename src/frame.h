#ifndef _FRAME_H
#define _FRAME_H

#include <stdio.h>
#include "actors.h"
#include "messages.h"
#include "intmap.h"
#include "pqueue.h"
#include "design.h"

struct frame_message {
	const struct message *msg;
	ssize_t interval;
};

struct frame_event {
	double time;
	ssize_t imsg;
};

/* dX[t,i] */
struct recv_frame {
	struct intmap jrecv_dxs;
};

struct frame {
	const struct design *design;
	double time;

	struct array observers;
	struct array frame_messages;
	ssize_t cur_fmsg;
	struct pqueue events;

	struct recv_frame *recv_frames;
	struct frame_var *vars;
};

struct frame_callbacks {
	void (*message_add) (void *udata, struct frame * f,
			     const struct message * msg);
	void (*message_advance) (void *udata, struct frame * f,
				 const struct message * msg, ssize_t intvl);
	void (*recv_update) (void *udata, struct frame * f, ssize_t isend,
			     ssize_t jrecv, ssize_t dyn_index, double delta);
	void (*send_update) (void *udata, struct frame * f, ssize_t isend,
			     ssize_t dyn_index, double dx);
	void (*clear) (void *udata, struct frame * f);
};

struct frame_observer {
	void *udata;
	struct frame_callbacks callbacks;
};

/* create/destroy/clear */
void frame_init(struct frame *f, const struct design *design);
void frame_deinit(struct frame *f);
void frame_clear(struct frame *f);

/* time */
double frame_time(const struct frame *f);	// current time
double frame_next_time(const struct frame *f);	// next change
void frame_advance(struct frame *f, double time);	// advance time

/* add a message  */
void frame_add(struct frame *f, const struct message *msg);

/* current covariates */
const struct vector *frame_recv_dx(const struct frame *f, ssize_t isend,
				   ssize_t jrecv);
void frame_recv_update(struct frame *f, ssize_t isend, ssize_t jrecv,
		       ssize_t dyn_index, double delta);

// struct vector frame_send_x(struct frame *f, ssize_t isend);
// void frame_send_update(const struct frame *f, ssize_t isend,
//                     ssize_t dyn_index, double delta);

/* observers */
void frame_add_observer(struct frame *f, void *udata,
			const struct frame_callbacks *callbacks);
void frame_remove_observer(struct frame *f, void *udata);

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

// void frame_send_mul(double alpha, enum trans_op trans,
//                  const struct frame *f,
//                  const struct vector *x, double beta, struct vector *y);
//void frame_send_muls(double alpha, enum trans_op trans,
//                   const struct frame *f,
//                   const struct svector *x, double beta, struct vector *y);

/* inline function definitions */
static inline const struct design *frame_design(const struct frame *f)
{
	assert(f);
	return f->design;
}

#endif /* _FRAME_H */
