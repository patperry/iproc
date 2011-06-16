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

typedef struct _iproc_sloglik iproc_sloglik;

struct _iproc_sloglik {
	struct refcount refcount;

	struct model *model;
	ssize_t isend;

	double f;
	struct vector *grad;
	bool grad_cached;

	ssize_t nsend;
	struct svector *nrecv;
	struct svector *dxobs;

	double gamma;
	struct svector *dp;
	struct svector *dxbar;
};

iproc_sloglik *iproc_sloglik_new(struct model *model, ssize_t isend);
iproc_sloglik *iproc_sloglik_ref(iproc_sloglik * sll);
void iproc_sloglik_unref(iproc_sloglik * sll);

void iproc_sloglik_insert(iproc_sloglik * sll,
			  const struct frame *f, ssize_t jrecv);
void iproc_sloglik_insertm(iproc_sloglik * sll,
			   const struct frame *f, ssize_t *jrecv, ssize_t n);

double iproc_sloglik_value(iproc_sloglik * sll);
struct vector *iproc_sloglik_grad(iproc_sloglik * sll);

#endif /* _IPROC_LOGLIK_H */
