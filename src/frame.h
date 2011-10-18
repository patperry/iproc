#ifndef _FRAME_H
#define _FRAME_H

#include <math.h>
#include <stdio.h>
#include "messages.h"
#include "pqueue.h"
#include "design.h"

struct frame_message {
	const struct message *message;
	size_t interval;
};

struct frame_event {
	double time;
	size_t imsg;
};

struct frame_actor {
	struct array message_ixs;
};

/* dX[t,i] */
struct recv_frame {
	size_t *active;	// active jrecv
	struct vector *dx;
	size_t nactive;
	size_t nactive_max;
};

struct frame {
	size_t nsend, nrecv;
	int has_loops;

	double *intvls;
	size_t nintvl;

	struct design recv_design;
	double time;

	struct array observers;
	struct array frame_messages;
	struct frame_actor *senders;
	struct frame_actor *receivers;
	size_t cur_fmsg;
	struct pqueue events;

	struct recv_frame *recv_frames;
};

struct frame_callbacks {
	void (*message_add) (void *udata, struct frame * f,
			     const struct message * msg);
	void (*message_advance) (void *udata, struct frame * f,
				 const struct message * msg, size_t intvl);
	void (*recv_update) (void *udata, struct frame * f, size_t isend,
			     size_t jrecv, const struct svector * delta);
	void (*send_update) (void *udata, struct frame * f, size_t isend,
			     size_t dyn_index, double dx);
	void (*clear) (void *udata, struct frame * f);
};

struct frame_observer {
	void *udata;
	struct frame_callbacks callbacks;
};

/* create/destroy/clear */
void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl);				
void frame_deinit(struct frame *f);
void frame_clear(struct frame *f);

/* properties */
static inline const double *frame_intervals(const struct frame *f);
static inline size_t frame_interval_count(const struct frame *f);

static inline struct design *frame_recv_design(const struct frame *f);

/* time */
static inline double frame_time(const struct frame *f);	// current time
static inline double frame_next_time(const struct frame *f);	// next change
void frame_advance(struct frame *f, double time);	// advance time

/* messages */
static inline size_t frame_messages_count(const struct frame *f);
static inline struct frame_message *frame_messages_item(const struct frame *f,
							size_t i);
void frame_add(struct frame *f, const struct message *msg);

/* actors */
static inline size_t frame_send_count(const struct frame *f);
static inline size_t frame_recv_count(const struct frame *f);
static inline int frame_has_loops(const struct frame *f);
static inline void frame_get_send_messages(const struct frame *f, size_t isend,
					   size_t **imsg, size_t *nmsg);
static inline void frame_get_recv_messages(const struct frame *f, size_t irecv,
					   size_t **imsg, size_t *nmsg);

/* current covariates */
void frame_recv_get_dx(const struct frame *f, size_t isend,
		       struct vector **dxp, size_t **activep,
		       size_t *nactivep);
const struct vector *frame_recv_dx(const struct frame *f, size_t isend,
				   size_t jrecv);
void frame_recv_update(struct frame *f, size_t isend, size_t jrecv,
		       const struct svector *delta);

// struct vector frame_send_x(struct frame *f, size_t isend);
// void frame_send_update(const struct frame *f, size_t isend,
//                     size_t dyn_index, double delta);

/* observers */
void frame_add_observer(struct frame *f, void *udata,
			const struct frame_callbacks *callbacks);
void frame_remove_observer(struct frame *f, void *udata);

void frame_recv_mul(double alpha, enum blas_trans trans,
		    const struct frame *f, size_t isend,
		    const struct vector *x, double beta, struct vector *y);
void frame_recv_muls(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const struct svector *x, double beta, struct vector *y);

void frame_recv_dmul(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const struct vector *x, double beta, struct svector *y);
void frame_recv_dmuls(double alpha, enum blas_trans trans,
		      const struct frame *f, size_t isend,
		      const struct svector *x, double beta, struct vector *y);

// void frame_send_mul(double alpha, enum blas_trans trans,
//                  const struct frame *f,
//                  const struct vector *x, double beta, struct vector *y);
//void frame_send_muls(double alpha, enum blas_trans trans,
//                   const struct frame *f,
//                   const struct svector *x, double beta, struct vector *y);

/* inline function definitions */
const double *frame_intervals(const struct frame *f)
{
	return f->intvls;
}

size_t frame_interval_count(const struct frame *f)
{
	return f->nintvl;
}

struct design *frame_recv_design(const struct frame *f)
{
	assert(f);
	return &((struct frame *)f)->recv_design;
}

double frame_time(const struct frame *f)
{
	assert(f);
	return f->time;
}

double frame_next_time(const struct frame *f)
{
	assert(f);

	if (pqueue_count(&f->events)) {
		const struct frame_event *e = pqueue_top(&f->events);
		return e->time;
	} else {
		return INFINITY;
	}
}

size_t frame_messages_count(const struct frame *f)
{
	assert(f);
	return array_count(&f->frame_messages);
}

struct frame_message *frame_messages_item(const struct frame *f, size_t imsg)
{
	assert(f);
	assert(imsg < frame_messages_count(f));
	struct frame_message *base = array_to_ptr(&f->frame_messages);
	return &base[imsg];
}

size_t frame_send_count(const struct frame *f)
{
	assert(f);
	return f->nsend;
}

size_t frame_recv_count(const struct frame *f)
{
	assert(f);
	return f->nrecv;
}

int frame_has_loops(const struct frame *f)
{
	return f->has_loops;
}

void frame_get_send_messages(const struct frame *f, size_t isend,
			     size_t **imsg, size_t *nmsg)
{
	assert(f);
	assert(isend < frame_send_count(f));
	assert(imsg);
	assert(nmsg);

	const struct frame_actor *a = &f->senders[isend];
	*imsg = array_to_ptr(&a->message_ixs);
	*nmsg = array_count(&a->message_ixs);
}

void frame_get_recv_messages(const struct frame *f, size_t irecv,
			     size_t **imsg, size_t *nmsg)
{
	assert(f);
	assert(irecv < frame_recv_count(f));
	assert(imsg);
	assert(nmsg);

	const struct frame_actor *a = &f->receivers[irecv];
	*imsg = array_to_ptr(&a->message_ixs);
	*nmsg = array_count(&a->message_ixs);
}

#endif /* _FRAME_H */
