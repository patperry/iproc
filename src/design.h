#ifndef _DESIGN_H
#define _DESIGN_H

#include "actors.h"
#include "history.h"
#include "list.h"
#include "dyad-queue.h"
#include "matrix.h"
#include "messages.h"
#include "refcount.h"
#include "svector.h"
#include "vector.h"

/* I. Introduction
 * ---------------
 *
 * A "design" holds the covariates for the interaction process.  It can
 * be thought of as a J-by-p matrix, indexed by time and sender.
 * 
 *         X[t,i] = [ x[t,i,1], ..., x[t,i,J] ]^T;
 * 
 * that is, the j-th row of X[t,i] is x[t,i,j].
 *
 * We represent the design as
 *
 *         X[t,i] = X[0,i] + dX[t,i].
 *
 *
 * II. Static part
 * ---------------
 * The matrix X[0,i] is referred to as the static design matrix.  This matrix
 * has a special structure which enables efficent multiplication.  Let
 * s[i] be the vector of sender i's traits and let
 * R = [ r[1], ..., r[J] ]^T be the matrix of the receivers' traits.  Then,
 *
 *         X[0,i] = [ kronecker(R, s[i]^T)  0 ]
 *
 * Suppose that s is p-dimensional and R is J-by-q.  If Y is a p-by-q
 * matrix, then
 *
 *         kronecker(R, s[i]^T) * vec(Y) = R * (Y^T * s[i]).
 *
 * If y is J-dimensional, then
 *
 *         (kronecker(R, s[i]^T))^T * y = vec(s[i]^T * (R^T * y)^T).
 *
 *
 * III. Dynamic part
 * -----------------
 * The dynamic part dX[t,i] is computed by client-supplied functions.  These
 * functions in general will depend on the history of the process.  The only
 * restriction is that dX[0,i] = 0.  Write
 *
 *         dX[t,i] = [  0  dX{1}[t,i]  dX{2}[t,i]  ...  dX{K}[t,i] ]
 *
 * The presumption is that each dX{k}[t,i] is sparse, in the sence that
 * (a) most of the rows of dX{k}[t,i] are zero, and (b) the nonzero rows
 * of dX{k}[t,i] are sparse vectors.  For each dX{k}[t,i], the client
 * must provide the following:
 *
 *         A. `dim` : the number of columns in dX{k}[t,i];
 *
 *         B. `get_dxs(dX{k}, ctx[t,i], offset)` : a function that
 *            copies the rows of dX{k}[t,i] into
 *            ctx->dxs[:,offset:(offset+dim-1)];
 *
 *         C. `free(dX{k})` : a funciton to free the memory associated
 *            with dX{k}
 *
 *
 * In context ctx, the matrix dX[t,i] is represented as an array of
 * { j, dx[t,i,j] } pairs called ctx->dxs, where dx[t,i,j] is
 * a sparse vector.
 */

struct frame;			// forward declarations
struct design;
struct design_var;
struct frame_var;

struct var_type {
	uint8_t dyad_event_mask;
	
	bool (*init) (struct design_var *dv, const struct design *d);
	void (*deinit) (struct design_var *dv);
	
	bool (*frame_init) (struct frame_var *fv, struct frame *f);
	void (*frame_deinit) (struct frame_var *fv);
	void (*frame_clear) (struct frame_var *fv);
	
	bool (*handle_dyad) (struct frame_var *fv,
			     const struct dyad_event *e,
			     struct frame *f);
};

struct design_var {
	const struct var_type *type;
	ssize_t dim;
	ssize_t index;
	void *udata;
};

struct design {
	struct actors *senders;
	struct actors *receivers;
	struct vector intervals;
	ssize_t ireffects;
	ssize_t istatic, nstatic;
	ssize_t idynamic, ndynamic;
	ssize_t dim;
	bool reffects;
	bool loops;

	struct list vars;
	struct refcount refcount;
};

struct design *design_alloc(struct actors *senders, struct actors *receivers,
			    const struct vector *intervals);
struct design *design_ref(struct design *design);
void design_free(struct design *design);

bool design_init(struct design *design, struct actors *senders,
		 struct actors *receivers, const struct vector *intervals);
void design_deinit(struct design *design);

ssize_t design_dim(const struct design *design);
ssize_t design_nsender(const struct design *design);
ssize_t design_nreceiver(const struct design *design);
struct actors *design_senders(const struct design *design);
struct actors *design_receivers(const struct design *design);

const struct vector *design_intervals(const struct design *design);

bool design_add_var(struct design *design, const struct var_type *type);
void design_set_loops(struct design *design, bool loops);
bool design_loops(const struct design *design);
void design_set_reffects(struct design *design, bool reffects);
bool design_reffects(const struct design *design);

ssize_t design_traits_index(const struct design *design);
ssize_t design_reffects_index(const struct design *design);
ssize_t design_var_index(const struct design *design,
			 const struct var_type *type);

void design_mul0(double alpha,
		 enum trans_op trans,
		 const struct design *design,
		 ssize_t isend,
		 const struct vector *x, double beta, struct vector *y);
void design_muls0(double alpha,
		  enum trans_op trans,
		  const struct design *design,
		  ssize_t isend,
		  const struct svector *x, double beta, struct vector *y);

#endif /* _DESIGN_H */
