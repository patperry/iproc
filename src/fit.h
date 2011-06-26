#ifndef _RECV_FIT_H
#define _RECV_FIT_H

#include "bfgs.h"
#include "design.h"
#include "linalg.h"
#include "linesearch.h"
#include "matrix.h"
#include "messages.h"
#include "model.h"
#include "recv_loglik.h"
#include "vector.h"


struct recv_fit {
	const struct design *design;
	const struct messages *msgs;
	double penalty;

	struct frame frame;
	struct vector coefs;
	struct model model;
	struct recv_loglik loglik;
	
	/* optimization problem */
	struct vector scale; /* scaling for the coefiecnets */
	struct matrix ce_t;  /* equality constraints: ce * coef = be */
	struct vector be;    /* cont'd */

	/* optimization problem workspace */
	struct vector duals;  /* dual parameters */
	struct vector resid;  /* dual and primal residuals */
	struct matrix kkt;    /* KKT matrix */
	
	/* additional workspace */
	struct linesearch ls;
	struct linesearch_ctrl lsctrl;
	struct symeig eig;	
};


void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct vector *coefs0,
		   double penalty);
void recv_fit_deinit(struct recv_fit *fit);
bool recv_fit_step(struct recv_fit *fit);
bool recv_fit_converged(const struct recv_fit *fit);
const char *recv_fit_errmsg(const struct recv_fit *fit);


#endif /* _RECV_FIT_H */
