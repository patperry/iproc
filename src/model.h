#ifndef _IPROC_MODEL_H
#define _IPROC_MODEL_H

#include <stdbool.h>
#include "refcount.h"
#include "array.h"
#include "design.h"
#include "frame.h"
#include "intmap.h"
#include "vector.h"

/* I. Model
 * --------
 * For each (i,j) dyad, we have a vector of covariates at time t, denoted
 * x[t,i,j].  We represent the covariates as
 *
 *         x[t,i,j] = x[0,i,j] + dx[t,i,j];
 *
 * often, dx is sparse, specifically it is zero for most (t, i, j) triples.
 *
 * The basic model for the weight of dyad (i,j) at time t is
 *
 *         w[t,i,j] = exp{eta[t,i,j]} 1{j in J[t,i]},
 *
 * where
 *
 *         eta[t,i,j] = (beta, x[t,i,j]).
 *
 * Above, beta is the coeifficient vector, and J[t,i] is the receiver
 * set of sender i.  Currently, J[t,i] is only affected by whether
 * or not we include self-loops in the model.  By convention, we set
 * J[0,i] to include all recipients.
 *
 * At time t, the probability of j being the recipient, given that i is
 * sending a message is
 *
 *         p[t,i,j] = w[t,i,j] / W[t,i],
 *
 * where
 *
 *         W[t,i] = sum { w[t,i,j] : j in J[t,i] }.
 *
 *
 * II. Probabilities
 * -----------------
 *
 * We express the probability at time t by
 *
 *         p[t,i,j] = gamma * p[0,i,j] + dp[t,i,j]
 *
 * where
 *
 *         gamma = gamma[t,i] = W[0,i] / W[t,i].
 *
 * The benefit of this representaiton is that dp[t,i,j] is often
 * sparse.  To see this, note
 * 
 *         eta[t,i,j] = eta[0,i,j] + (beta, dx[t,i,j]),
 *           w[t,i,j] = w[0,i,j] * exp{(beta, dx[t,i,j])} 1{j in J[t,i]},
 *           p[t,i,j] = gamma * p[0,i,j] * exp{(beta, dx[t,i,j])} 1{j in J[t,i]}
 * 
 * Thus, if j is in J[t,i], then
 *
 *          dp[t,i,j] = (gamma * p[0,i,j]
 *                       * (exp{(beta, dx[t,i,j])} - 1));
 *
 * otherwise,
 *
 *          dp[t,i,j] = - gamma * p[0,i,j].
 *
 * For most (i,j) pairs, dx[t,i,j] = 0 and j is in J[t,i], so that this term
 * is zero.  We can simplify the expreesion for dp by defining
 * 
 *          deta[t,i,j] = (beta, dx[t,i,j]) if j is in J[t,i],
 *                        -Infinity         otherwise;
 *
 * then,
 *
 *          dp[t,i,j] = gamma * p[0,i,j] * (exp{deta[t,i,j]} - 1)
 * 
 * We express the log probability at time t by
 *
 *         log(p[t,i,j]) = log(gamma) + log(p[0,i,j]) + deta[t,i,j].
 *
 *
 * III. Covariate means
 * --------------------
 * The mean of the covariates for sender i at time t is defined as
 *
 *         xbar[t,i] = sum{ p[t,i,j] * x[t,i,j], j in J[t,i] };
 *
 * in matrix form
 *
 *         xbar[t,i] = X[t,i]^T * P[t,i]
 *                   = (X[0,i] + dX[t,i])^T * (gamma P[0,i] + dP[t,i])
 *                   = gamma * xbar[0,i]
 *                     + ( X[0,i])^T * dP[t,i]
 *                     + (dX[t,i])^T *  P[t,i]
 *
 * Defining dxbar[t,i] to be the latter term,
 *
 *         xbar[t,i] = gamma * xbar[0,i] + (X[0,i])^T * dP[t,i] + dxbar[t,i].
 */

typedef struct cohort_model iproc_group_model;
typedef struct _iproc_model iproc_model;
typedef struct _iproc_model_ctx iproc_model_ctx;

/* Two senders, i1 and i2, are in the same group if and only if their
 * covariates agree at time 0, i.e.
 *
 *   x[0,i1,j] = x[0,i2,j] for all j.
 *
 * Owing to the convention that J[0,i] includes all receivers,
 * all senders in the same group share the same values of
 * w[0,i,j], W[0,i,j], p[0,i,j], and xbar[0,i].
 */
struct cohort_model {
	double log_W0;
	struct vector p0;
	struct vector log_p0;
	struct vector xbar0;
};

struct _iproc_model {
	struct design *design;
	struct vector coefs;
	bool has_loops;
	struct intmap cohort_models;
	struct array ctxs;
	struct refcount refcount;
};

struct _iproc_model_ctx {
	iproc_model *model;
	const struct frame *frame;
	ssize_t isend;
	iproc_group_model *group;

	double gamma;
	double log_gamma;
	struct svector deta;
	struct svector dp;
	struct svector dxbar;

	struct refcount refcount;
};

iproc_model *iproc_model_new(struct design *design,
			     struct vector *coefs, bool has_loops);
iproc_model *iproc_model_ref(iproc_model * model);
void iproc_model_unref(iproc_model * model);

struct design *iproc_model_design(iproc_model * model);
struct vector *iproc_model_coefs(iproc_model * model);
bool iproc_model_has_loops(iproc_model * model);
ssize_t iproc_model_nsender(iproc_model * model);
ssize_t iproc_model_nreceiver(iproc_model * model);
ssize_t iproc_model_dim(iproc_model * model);

/* Initial probability, and expectations, without adjustment for self-loops. */
iproc_group_model *iproc_model_send_group(iproc_model * model, ssize_t isend);
struct vector *iproc_model_logprobs0(iproc_model * model, ssize_t isend);
struct vector *iproc_model_probs0(iproc_model * model, ssize_t isend);
struct vector *iproc_model_mean0(iproc_model * model, ssize_t isend);

iproc_model_ctx *iproc_model_ctx_new(iproc_model * model,
				     const struct frame *f, ssize_t isend);
iproc_model_ctx *iproc_model_ctx_ref(iproc_model_ctx * ctx);
void iproc_model_ctx_unref(iproc_model_ctx * ctx);

void iproc_model_ctx_set(iproc_model_ctx * ctx, const struct frame *f,
			 ssize_t isend);

ssize_t iproc_model_ctx_nreceiver(iproc_model_ctx * ctx);
double iproc_model_ctx_prob(iproc_model_ctx * ctx, ssize_t jrecv);
double iproc_model_ctx_logprob(iproc_model_ctx * ctx, ssize_t jrecv);
void iproc_model_ctx_get_probs(iproc_model_ctx * ctx, struct vector *probs);
void iproc_model_ctx_get_logprobs(iproc_model_ctx * ctx,
				  struct vector *logprobs);

/*
void                iproc_vector_acc_diffmean (struct vector    *dst_vector,
                                               double           scale,
                                               iproc_model_ctx *ctx);
*/

#endif /* _IPROC_MODEL_H */
