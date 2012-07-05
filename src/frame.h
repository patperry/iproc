#ifndef _FRAME_H
#define _FRAME_H

#include <inttypes.h> // imaxdiv, imaxdiv_t
#include <math.h>
#include <stdio.h>
#include "sblas.h"
#include "messages.h"
#include "pqueue.h"
#include "history.h"
#include "design.h"

/* dX[t,i] */
struct recv_frame {
	struct vpattern active;
};

struct frame {
	struct history history;
	size_t nsend, nrecv;
	struct design send_design;
	struct design recv_design;
	struct design dyad_design;

	int has_loops;

	double *intvls;
	size_t nintvl;

	double time;

	struct frame_observer *observers;
	size_t nobs, nobs_max;
	struct recv_frame *recv_frames;
};

struct frame_callbacks {
	void (*message_add) (void *udata, struct frame * f,
			     const struct message * msg);
	void (*message_advance) (void *udata, struct frame * f,
				 const struct message * msg, size_t intvl);
	void (*recv_update) (void *udata, struct frame * f, size_t isend,
			     size_t jrecv, const double *delta,
			     const struct vpattern *pat);
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

static inline struct history *frame_history(const struct frame *f);
static inline struct design *frame_send_design(const struct frame *f);
static inline struct design *frame_recv_design(const struct frame *f);
static inline struct design *frame_dyad_design(const struct frame *f);

/* time */
static inline double frame_time(const struct frame *f);	// current time
static inline double frame_next_time(const struct frame *f);	// next change
void frame_advance(struct frame *f, double time);	// advance time

/* messages */
void frame_add(struct frame *f, const struct message *msg);

/* actors */
static inline size_t frame_send_count(const struct frame *f);
static inline size_t frame_recv_count(const struct frame *f);
static inline size_t frame_dyad_ix(const struct frame *f, size_t isend, size_t jrecv);
static inline void frame_get_dyad(const struct frame *f, size_t ix, size_t *pisend, size_t *pjrecv);
static inline int frame_has_loops(const struct frame *f);


/* current covariates */
void frame_recv_get_dx(const struct frame *f, size_t isend,
		       const double **dxp, const size_t **activep,
		       size_t *nactivep);
const double *frame_recv_dx(const struct frame *f, size_t isend, size_t jrecv);
//void frame_recv_update(struct frame *f, size_t isend, size_t jrecv,
//		       const double *delta, const struct vpattern *pat);

// const double *frame_send_x(struct frame *f, size_t isend);
// void frame_send_update(const struct frame *f, size_t isend,
//                     size_t dyn_index, double delta);

/* observers */
void frame_add_observer(struct frame *f, void *udata,
			const struct frame_callbacks *callbacks);
void frame_remove_observer(struct frame *f, void *udata);

/*
void frame_recv_mul(double alpha, enum blas_trans trans,
		    const struct frame *f, size_t isend,
		    const double *x, double beta, double *y);
void frame_recv_muls(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, const struct vpattern *pat,
		     double beta, double *y);

void frame_recv_dmul(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, double beta, double *y);
void frame_recv_dmuls(double alpha, enum blas_trans trans,
		      const struct frame *f, size_t isend,
		      const double *x, const struct vpattern *pat, double beta,
		      double *y);
*/

// void frame_send_mul(double alpha, enum blas_trans trans,
//                  const struct frame *f,
//                  const double *x, double beta, double *y);
//void frame_send_muls(double alpha, enum blas_trans trans,
//                   const struct frame *f,
//                   const double *x, const struct vpattern *pat,
//                   double beta, struct double *y);

/* inline function definitions */
const double *frame_intervals(const struct frame *f)
{
	return history_intervals(&f->history);
}

size_t frame_interval_count(const struct frame *f)
{
	return history_interval_count(&f->history);
}

struct history *frame_history(const struct frame *f)
{
	return &((struct frame *)f)->history;
}

struct design *frame_send_design(const struct frame *f)
{
	assert(f);
	return &((struct frame *)f)->send_design;
}

struct design *frame_recv_design(const struct frame *f)
{
	assert(f);
	return &((struct frame *)f)->recv_design;
}

struct design *frame_dyad_design(const struct frame *f)
{
	assert(f);
	return &((struct frame *)f)->dyad_design;
}

double frame_time(const struct frame *f)
{
	assert(f);
	return history_time(&f->history);
}

double frame_next_time(const struct frame *f)
{
	assert(f);
	return history_next_time(&f->history);
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

size_t frame_dyad_ix(const struct frame *f, size_t isend, size_t jrecv)
{
	assert(0 <= isend && isend < frame_send_count(f));
	assert(0 <= jrecv && jrecv < frame_recv_count(f));

	size_t nrecv = frame_recv_count(f);
	size_t ix = jrecv + isend * nrecv;
	return ix;
}

void frame_get_dyad(const struct frame *f, size_t ix, size_t *pisend, size_t *pjrecv)
{
	size_t nrecv = frame_recv_count(f);
	imaxdiv_t res = imaxdiv(ix, nrecv);
	size_t isend = res.quot;
	size_t jrecv = res.rem;
	
	*pisend = isend;
	*pjrecv = jrecv;
}


int frame_has_loops(const struct frame *f)
{
	return f->has_loops;
}

#endif /* _FRAME_H */
