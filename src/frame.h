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
#include "design2.h"

struct dyad {
	size_t isend;
	size_t jrecv;
};


struct frame {
	struct history history;
	size_t nsend, nrecv;
	struct design send_design;
	struct design recv_design;
	struct design2 dyad_design;

	int has_loops;
};


/* create/destroy/clear */
void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl);
void frame_deinit(struct frame *f);

/* properties */
static inline struct history *frame_history(const struct frame *f);
static inline struct design *frame_send_design(const struct frame *f);
static inline struct design *frame_recv_design(const struct frame *f);
static inline struct design2 *frame_dyad_design(const struct frame *f);

/* actors */
static inline size_t frame_send_count(const struct frame *f);
static inline size_t frame_recv_count(const struct frame *f);
static inline size_t frame_dyad_ix(const struct frame *f, size_t isend, size_t jrecv);
static inline struct dyad frame_ix_dyad(const struct frame *f, size_t ix);
static inline int frame_has_loops(const struct frame *f);




/* inline function definitions */
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

struct design2 *frame_dyad_design(const struct frame *f)
{
	assert(f);
	return &((struct frame *)f)->dyad_design;
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
	assert(isend < frame_send_count(f));
	assert(jrecv < frame_recv_count(f));

	size_t nrecv = frame_recv_count(f);
	size_t ix = jrecv + isend * nrecv;
	return ix;
}

struct dyad frame_ix_dyad(const struct frame *f, size_t ix)
{
	size_t nrecv = frame_recv_count(f);
	imaxdiv_t res = imaxdiv(ix, nrecv);

	struct dyad dyad;
	dyad.isend = (size_t)res.quot;
	dyad.jrecv = (size_t)res.rem;
	
	return dyad;
}

int frame_has_loops(const struct frame *f)
{
	return f->has_loops;
}

#endif /* _FRAME_H */
