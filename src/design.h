#ifndef _DESIGN_H
#define _DESIGN_H

#include "actors.h"
#include "history.h"
#include "darray.h"
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

struct frame; // forward declaration
struct dyad_var;

struct dyad_var {
	ssize_t dim;
	uint8_t dyad_event_mask;
	bool (*handle_dyad_event) (struct dyad_var *v,
				   const struct dyad_event *e,
				   struct frame *f,
				   ssize_t index);
	bool (*get_jrecv_dxs) (struct dyad_var *v,
			       struct frame *f,
			       ssize_t index);
};


struct design_dyad_var {
	ssize_t index;
	struct dyad_var *var;
};

struct design {
	struct actors *senders;
	struct actors *receivers;
	const struct vector *intervals;
	ssize_t ireffects, nreffects;
	ssize_t istatic, nstatic;
	ssize_t idynamic, ndynamic;
	ssize_t dim;
	bool has_reffects;
	bool has_loops;
	
	struct darray design_dyad_vars;

	/* deprecated */
	struct darray vars;
	struct darray ctxs;
	struct darray svectors;
	struct refcount refcount;
};

struct design *design_alloc(struct actors *senders, struct actors *receivers,
			    bool has_reffects, bool has_loops,
			    const struct vector *intervals);
struct design *design_ref(struct design * design);
void design_free(struct design * design);

bool design_init(struct design *design, struct actors *senders,
		 struct actors *receivers, bool has_reffects, bool has_loops,
		 const struct vector *intervals);
void design_deinit(struct design *design);


ssize_t design_dim(const struct design * design);
ssize_t design_nsender(const struct design * design);
ssize_t design_nreceiver(const struct design * design);
struct actors *design_senders(const struct design * design);
struct actors *design_receivers(const struct design * design);
bool design_has_reffects(const struct design *design);
bool design_has_loops(const struct design *design);
const struct vector *design_intervals(const struct design *design);


void design_mul0(double alpha,
		       enum trans_op trans,
		       const struct design * design,
		       ssize_t isend,
		       const struct vector *x, double beta, struct vector *y);
void design_muls0(double alpha,
			enum trans_op trans,
			const struct design * design,
			ssize_t isend,
			const struct svector *x, double beta, struct vector *y);


ssize_t design_add_dyad_var(struct design *design, struct dyad_var *var);



/********** DEPRECATED **********/

typedef struct _iproc_design_var iproc_design_var;
typedef struct _iproc_design_ctx iproc_design_ctx;
typedef struct _iproc_design_dx iproc_design_dx;	// private


/* dX[t,i] */
struct _iproc_design_ctx {
	struct design *design;
	struct history *history;
	ssize_t isend;
	struct darray dxs;
	struct refcount refcount;
};

/* dX{k}[t,i] */
struct _iproc_design_var {
	ssize_t dim;
	void (*get_dxs) (iproc_design_var * var, iproc_design_ctx * ctx,
			 ssize_t offset);
	void (*free) (iproc_design_var * var);
};

/* dx[t,i,j] */
struct _iproc_design_dx {
	ssize_t jrecv;		// NOTE: layout is important here
	struct svector *dx;
};

// NOTE: a call to `append` invalidates existing design_ctxs
void iproc_design_append(struct design * design, iproc_design_var * var);


void iproc_design_var_init(iproc_design_var * var,
			   ssize_t dim,
			   void (*get_dxs) (iproc_design_var *,
					    iproc_design_ctx * ctx,
					    ssize_t),
			   void (*free) (iproc_design_var *));

iproc_design_ctx *iproc_design_ctx_new(struct design * design,
				       ssize_t isend, struct history * h);
iproc_design_ctx *iproc_design_ctx_ref(iproc_design_ctx * ctx);
void iproc_design_ctx_unref(iproc_design_ctx * ctx);

void iproc_design_ctx_mul(double alpha,
			  enum trans_op trans,
			  iproc_design_ctx * ctx,
			  struct vector *x, double beta, struct vector *y);
void iproc_design_ctx_muls(double alpha,
			   enum trans_op trans,
			   iproc_design_ctx * ctx,
			   struct svector *x, double beta, struct vector *y);

struct svector *iproc_design_ctx_dx(iproc_design_ctx * ctx,
				    ssize_t jrecv, bool null_ok);

/* number of nonzero rows in dX */
ssize_t iproc_design_ctx_nnz(iproc_design_ctx * ctx);

/* inz-th nonzero row of dX */
struct svector *iproc_design_ctx_nz(iproc_design_ctx * ctx,
				    ssize_t inz, ssize_t *jrecv);

void iproc_design_ctx_dmul(double alpha,
			   enum trans_op trans,
			   const iproc_design_ctx * ctx,
			   const struct vector *x, double beta,
			   struct svector *y);
void iproc_design_ctx_dmuls(double alpha, enum trans_op trans,
			    const iproc_design_ctx * ctx,
			    const struct svector *x, double beta,
			    struct svector *y);



#endif /* _DESIGN_H */
