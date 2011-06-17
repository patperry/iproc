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

struct recv_sloglik {
	struct model *model;
	ssize_t isend;

	double f;
	struct vector grad;
	bool grad_cached;

	ssize_t nsend;
	struct svector nrecv;
	struct vector dxobs;

	double gamma;
	struct svector dp;
	struct vector dxbar;
};

void recv_sloglik_init(struct recv_sloglik *ll, const struct model *model, ssize_t isend);
void recv_sloglik_deinit(struct recv_sloglik *ll);

struct recv_sloglik *recv_sloglik_alloc(const struct model *model, ssize_t isend);
void recv_sloglik_free(struct recv_sloglik *ll);

void recv_sloglik_add(struct recv_sloglik *ll,
		     const struct frame *f, ssize_t *jrecv, ssize_t n);

double recv_sloglik_value(const struct recv_sloglik *ll);
void recv_sloglik_axpy_grad(double alpha, const struct recv_sloglik *ll, struct vector *y);

#endif /* _IPROC_SLOGLIK_H */
