#ifndef _IPROC_SLOGLIK_H
#define _IPROC_SLOGLIK_H

#include <stdbool.h>

#include "array.h"
#include "frame.h"
#include "history.h"
#include "model.h"
#include "refcount.h"
#include "svector.h"

/* I. Value
 * --------
 * The value of the negative log partial likelihood for sender i is
 *
 *         f = -sum{ log(p[t,i,j]) }
 *
 * where the sum is over all observed message triples (t,i,j).
 * 
 *
 * II. Gradient
 * ------------
 * The gradient of the negative log partial likelihood is
 *
 *         g = sum{ xbar[t,i] - x[t,i,j] }.
 *
 * We use the representation
 *
 *         xbar[t,i] = gamma * xbar[0,i]
 *                     + ( X[0,i])^T * dP[t,i]
 *                     + dxbar[t,i],
 *
 * so that
 *
 *         g = [ sum{gamma[t,i]} * xbar[0,i]
 *               + ( X[0,i])^T * sum{dP[t,i]}
 *               + sum{dxbar[t,i]} ]
 *             -
 *             [ (X[0,i])^T n[i] + sum{dx[t,i,j]} ],
 *
 * where n[i] is a vector whose jth entry is the number of messages
 * sent from i to j.
 *
 *
 * III. Hessian
 * ------------
 * The Hessian of the negative log partial likelihood is
 *
 *         H = sum{ ... }
 *
 */

struct recv_sloglik_score {
	const struct vector *mean0;
	struct svector nrecv;
	struct vector mean_obs_dx;
	double gamma;
	struct vector dp;	// p - gamma * p0
	struct vector mean_dx;	// dx' * p
};

struct recv_sloglik_imat {
	const struct matrix *imat0;
	double gamma2;		// gamma * (1 - gamma)
	struct vector gamma_dp; // gamma * dp
	struct vector gamma_mean_dx; // gamma * dx' * p
	struct matrix dx_p;     // dx' * diag(p)
	struct matrix mean_dx_dp; // (dx' * p) * dp'
	struct matrix dp2;	// diag(dp) - dp * dp'
	struct matrix var_dx;	// dx' * [diag(p) - p p'] * dx
};

struct recv_sloglik {
	struct model *model;
	ssize_t isend;

	ssize_t n_last, n;
	double dev_last, dev_avg;	
	struct array active;	
	struct recv_sloglik_score score_last, score_avg;
	struct recv_sloglik_imat imat_last, imat_avg;		
};

void recv_sloglik_init(struct recv_sloglik *ll, const struct model *model, ssize_t isend);
void recv_sloglik_deinit(struct recv_sloglik *ll);

void recv_sloglik_add(struct recv_sloglik *ll,
		      const struct frame *f, const ssize_t *jrecv, ssize_t n);
void recv_sloglik_clear(struct recv_sloglik *ll);

ssize_t recv_sloglik_count(const struct recv_sloglik *sll);
double recv_sloglik_avg_dev(const struct recv_sloglik *sll);
void recv_sloglik_axpy_avg_mean(double alpha, const struct recv_sloglik *sll, struct vector *y);
void recv_sloglik_axpy_avg_score(double alpha, const struct recv_sloglik *sll, struct vector *y);
void recv_sloglik_axpy_avg_imat(double alpha, const struct recv_sloglik *sll, struct matrix *y);

ssize_t recv_sloglik_last_count(const struct recv_sloglik *sll);
double recv_sloglik_last_dev(const struct recv_sloglik *sll);
void recv_sloglik_axpy_last_mean(double alpha, const struct recv_sloglik *sll, struct vector *y);
void recv_sloglik_axpy_last_score(double alpha, const struct recv_sloglik *sll, struct vector *y);
void recv_sloglik_axpy_last_imat(double alpha, const struct recv_sloglik *sll, struct matrix *y);


#endif /* _IPROC_SLOGLIK_H */
