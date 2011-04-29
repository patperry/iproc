#ifndef _FRAME_H
#define _FRAME_H

#include "actors.h"
#include "intmap.h"
#include "pqueue.h"

struct frame {
	struct actors *senders;
	struct actors *receivers;
	bool has_reffects;
	ssize_t ireffects, nreffects;
	ssize_t istatic, nstatic;
	ssize_t idynamic, ndynamic;
	ssize_t dim;
	
	double time;
	struct darray dyad_vars;
	struct pqueue dyad_var_diffs;
	struct intmap isend_frames; // (j, dX[t,i) pairs; dX is a 'struct send_frame'
};

/* dX[t,i] */
struct send_frame {
	ssize_t isend;
	struct intmap jrecv_dxs; // (j, dx[t,i,j]) pairs; dx is a 'struct svector'
};

struct message {
	ssize_t from;
	ssize_t *to;
	ssize_t nto;
	intptr_t attr;
};

struct dyad_var {
	ssize_t dim;
	ssize_t offset;
	bool (*insert) (struct dyad_var *v, const struct message *msg, struct frame *f);
	void (*deinit) (struct dyad_var *v);
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

void design_var_init(struct design_var *v, ssize_t dim,
		     bool (*insert) (struct dyad_var *v, const struct message *msg, struct frame *f),
		     void (*deinit) (struct dyad_var *v));
void design_var_deinit(struct design_var *v);

/* create/destroy */
bool frame_init(struct frame *f, const struct actors *senders,
		const struct actors *receivers);
void frame_deinit(struct frame *f);

/* append a new variable to the design matrix */
bool frame_add_var(struct frame *f, struct design_var *x);

bool frame_add_dyad_event(struct frame *f, double t,
			  const struct design_dyad_var_diff delta);
bool frame_add_send_event(struct frame *f, double t,
			  const struct design_send_var_diff delta);


/* record a message event */
bool frame_insert(struct frame *f, const struct message *msg);

/* advance time */
void frame_advance_to(const struct frame *f, double t, struct frame_diff *diff);

/* time of the next change */
double frame_time(const struct frame *f);
double frame_next_update(const struct frame *f);



void frame_mul0(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct vector *x, double beta, struct vector *y);
void frame_muls0(double alpha, enum trans_op trans,
		 const struct frame *f, ssize_t isend,
		 const struct svector *x, double beta, struct vector *y);

void frame_mul(double alpha, enum trans_op trans,
	       const struct frame *f, ssize_t isend,
	       const struct vector *x, double beta, struct vector *y);
void frame_muls(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct svector *x, double beta, struct vector *y);

void frame_dmul(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct vector *x, double beta, struct svector *y);
void frame_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f, ssize_t isend,
		 const struct svector *x, double beta, struct svector *y);



#endif /* _FRAME_H */