#ifndef FRAME_H
#define FRAME_H

#include "history.h"
#include "design.h"
#include "design2.h"
#include <stddef.h> /* size_t */


struct frame {
	struct history history;
	struct design send_design;
	struct design recv_design;
	struct design2 dyad_design;
	int has_loops;
};


/* create/destroy */
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
	return design_count(&f->send_design);
}

size_t frame_recv_count(const struct frame *f)
{
	return design_count(&f->recv_design);
}

int frame_has_loops(const struct frame *f)
{
	return f->has_loops;
}

#endif /* FRAME_H */
