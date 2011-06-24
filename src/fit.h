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
	struct model model;
	struct recv_loglik loglik;
	
	double f;
	struct vector grad;
	struct matrix imat_evec;
	struct vector imat_eval;
	ssize_t imat_rank;
	struct symeig eig;
	
	struct linesearch ls;
	struct linesearch_ctrl lsctrl;
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
