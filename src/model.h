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
	ssize_t isend0;		// canonical sender
	double max_eta0;
	double log_W0;		// log_W0 - max_eta0
	struct vector eta0;
	struct vector p0;
	struct vector mean0;
	struct matrix imat0;

	/* debug */
#ifndef NDEBUG
	double W0;
	struct vector w0;
#endif
};

struct recv_model {
	struct model *model;
	ssize_t isend;
	struct cohort_model *cohort;
	double gamma;
	double log_W;		// log_W - scale
	double scale;
	struct svector deta;
	struct array active;
};

struct model {
	struct frame *frame;
	struct vector recv_coefs;
	struct intmap cohort_models;
	struct intmap recv_models;
	struct refcount refcount;
};

void model_init(struct model *model,
		struct frame *f, const struct vector *recv_coefs);
void model_deinit(struct model *model);
struct model *model_alloc(struct frame *f, const struct vector *recv_coefs);
struct model *model_ref(struct model *model);
void model_free(struct model *model);

const struct frame *model_frame(const struct model *model);
const struct design *model_design(const struct model *model);
const struct vector *model_recv_coefs(const struct model *model);
ssize_t model_send_count(const struct model *model);
ssize_t model_recv_count(const struct model *model);
ssize_t model_recv_dim(const struct model *model);

void model_set_recv_coefs(struct model *m, const struct vector *recv_coefs);

struct recv_model *model_recv_model(const struct model *m, ssize_t isend);
ssize_t recv_model_count(const struct recv_model *rm);
ssize_t recv_model_dim(const struct recv_model *rm);

/* Initial probability, and expectations, without adjustment for self-loops. */
static inline double recv_model_prob0(const struct recv_model *rm,
				      ssize_t jrecv);

struct vector *recv_model_logweight0(const struct recv_model *rm);
struct vector *recv_model_probs0(const struct recv_model *rm);
struct vector *recv_model_mean0(const struct recv_model *rm);
struct matrix *recv_model_imat0(const struct recv_model *rm);

double recv_model_logprob(const struct recv_model *rm, ssize_t jrecv);
double recv_model_prob(const struct recv_model *rm, ssize_t jrecv);
void recv_model_axpy_probs(double alpha, const struct recv_model *rm,
			   struct vector *y);

static inline double recv_model_prob0(const struct recv_model *rm,
				      ssize_t jrecv)
{
	assert(rm);
	assert(0 <= jrecv && jrecv < recv_model_count(rm));

	const struct vector *p0 = &rm->cohort->p0;
	return vector_item(p0, jrecv);
}

#endif /* _IPROC_MODEL_H */
