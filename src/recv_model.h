#ifndef _IPROC_MODEL_H
#define _IPROC_MODEL_H

#include <stddef.h>
#include "design.h"
#include "frame.h"
#include "mlogit.h"
#include "mlogitaug.h"

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

struct recv_coefs {
	double *all;
	size_t dim;
	struct coefs recv;
	struct coefs2 dyad;
};

struct recv_model_cohort {
	struct mlogit mlogit;
	double *eta0;
};

struct recv_model_sender {
	struct mlogitaug mlogitaug;
};

struct recv_model {
	struct frame *frame;
	struct recv_coefs coefs;
	struct recv_model_cohort *cohort_models;
	struct recv_model_sender *sender_models;
};


void recv_coefs_init(struct recv_coefs *c, const struct frame *f);
void recv_coefs_deinit(struct recv_coefs *c);

void recv_model_init(struct recv_model *model,
		     struct frame *f,
		     const struct recv_coefs *coefs);
void recv_model_deinit(struct recv_model *model);

struct frame *recv_model_frame(const struct recv_model *model);
const struct design *recv_model_design(const struct recv_model *model);
const size_t *recv_model_cohorts(const struct recv_model *model);
const struct recv_coefs *recv_model_coefs(const struct recv_model *model);
size_t recv_model_send_count(const struct recv_model *model);
size_t recv_model_cohort_count(const struct recv_model *model);
size_t recv_model_cohort(const struct recv_model *model, size_t isend);
size_t recv_model_count(const struct recv_model *model);
size_t recv_model_dim(const struct recv_model *model);

void recv_model_set_coefs(struct recv_model *m, const struct recv_coefs *coefs);

double recv_model_lprob(const struct recv_model *m, size_t isend,
			  size_t jrecv);
double recv_model_prob(const struct recv_model *m, size_t isend,
		       size_t jrecv);
void recv_model_axpy_probs(double alpha, const struct recv_model *m,
			   size_t isend, double *y);

double recv_model_psi(const struct recv_model *m, size_t isend);


#if 0

/* Initial probability, and expectations, without adjustment for self-loops. */
double recv_model_logsumwt0(const struct recv_model *m, size_t c);
double *recv_model_logwts0(const struct recv_model *m, size_t c);
double *recv_model_probs0(const struct recv_model *m, size_t c);

double recv_model_prob0(const struct recv_model *m, size_t c, size_t jrecv);
double *recv_model_mean0(const struct recv_model *m, size_t c);
double *recv_model_imat0(const struct recv_model *m, size_t c);




void recv_model_get_active(const struct recv_model *m, size_t isend,
			   size_t **jrecv, size_t *n);
double recv_model_invgrow(const struct recv_model *m, size_t isend);

#endif

#endif /* _IPROC_MODEL_H */
