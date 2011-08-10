#ifndef _DESIGN_H
#define _DESIGN_H

#include "actors.h"
#include "history.h"
#include "array.h"
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

struct var_type; // forward declaration

struct design {
	struct actors *senders;
	struct actors *receivers;
	const struct matrix *traits;
	bool loops;

	struct vector intervals;

	struct array send_vars;
	bool seffects;
	ssize_t iseffects;
	ssize_t isstatic, nsstatic;
	ssize_t isdynamic, nsdynamic;
	ssize_t sdim;

	struct array recv_vars;
	bool reffects;
	ssize_t ireffects;
	ssize_t irstatic, nrstatic;
	ssize_t irdynamic, nrdynamic;
	ssize_t rdim;

	struct refcount refcount;
};

struct design *design_alloc(struct actors *senders, struct actors *receivers,
			    const struct matrix *traits,
			    const struct vector *intervals);
struct design *design_ref(struct design *design);
void design_free(struct design *design);

void design_init(struct design *design, struct actors *senders,
		 struct actors *receivers, const struct matrix *traits,
		 const struct vector *intervals);
void design_deinit(struct design *design);

static inline ssize_t design_send_dim(const struct design *design);
static inline ssize_t design_recv_dim(const struct design *design);
static inline ssize_t design_send_count(const struct design *design);
static inline ssize_t design_recv_count(const struct design *design);
static inline struct actors *design_senders(const struct design *design);
static inline struct actors *design_receivers(const struct design *design);
static inline const struct matrix *design_traits(const struct design *design);
static inline const struct vector *design_intervals(const struct design
						    *design);

static inline bool design_loops(const struct design *design);
void design_set_loops(struct design *design, bool loops);

bool design_send_effects(const struct design *design);
void design_set_send_effects(struct design *design, bool seffects);
void design_add_send_var(struct design *design, const struct var_type *type);
ssize_t design_send_traits_index(const struct design *design);
ssize_t design_send_effects_index(const struct design *design);
ssize_t design_send_var_index(const struct design *design,
			      const struct var_type *type);

/*
void design_send_mul0(double alpha,
		      enum trans_op trans,
		      const struct design *design,
		      const struct vector *x, double beta, struct vector *y);
void design_send_muls0(double alpha,
		       enum trans_op trans,
		       const struct design *design,
		       const struct svector *x, double beta, struct vector *y);
*/

bool design_recv_effects(const struct design *design);
void design_set_recv_effects(struct design *design, bool reffects);
void design_add_recv_var(struct design *design, const struct var_type *type);
ssize_t design_recv_effects_index(const struct design *design);
ssize_t design_recv_traits_index(const struct design *design);
ssize_t design_recv_var_index(const struct design *design,
			      const struct var_type *type);
static inline ssize_t design_recv_dyn_index(const struct design *design);
static inline ssize_t design_recv_dyn_dim(const struct design *design);

void design_recv_mul0(double alpha,
		      enum trans_op trans,
		      const struct design *design,
		      const struct vector *x, double beta, struct vector *y);
void design_recv_muls0(double alpha,
		       enum trans_op trans,
		       const struct design *design,
		       const struct svector *x, double beta, struct vector *y);

/* inline funciton definitions */
ssize_t design_send_dim(const struct design *design)
{
	assert(design);
	return design->sdim;
}

ssize_t design_recv_dim(const struct design *design)
{
	assert(design);
	return design->rdim;
}

ssize_t design_send_count(const struct design *design)
{
	assert(design);
	const struct actors *senders = design_senders(design);
	return actors_count(senders);
}

ssize_t design_recv_count(const struct design *design)
{
	assert(design);
	const struct actors *receivers = design_receivers(design);
	return actors_count(receivers);
}

struct actors *design_senders(const struct design *design)
{
	assert(design);
	return design->senders;
}

struct actors *design_receivers(const struct design *design)
{
	assert(design);
	return design->receivers;
}

const struct matrix *design_traits(const struct design *design)
{
	assert(design);
	return design->traits;
}

const struct vector *design_intervals(const struct design *design)
{
	assert(design);
	return &design->intervals;
}

bool design_loops(const struct design *design)
{
	assert(design);
	return design->loops;
}

ssize_t design_recv_dyn_index(const struct design *design)
{
	assert(design);
	return design->irdynamic;
}

ssize_t design_recv_dyn_dim(const struct design *design)
{
	assert(design);
	return design->nrdynamic;
}

#endif /* _DESIGN_H */
