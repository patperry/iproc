#ifndef _IPROC_FIT_H
#define _IPROC_FIT_H

#include <stdbool.h>
#include <float.h>
#include <math.h>

#include "bfgs.h"
#include "design.h"
#include "recv_loglik.h"
#include "matrix.h"
#include "messages.h"
#include "model.h"
#include "vector.h"


struct recv_fit {
	const struct design *design;
	const struct messages *msgs;
	double penalty;

	struct frame frame;	
	struct model model;
	struct recv_loglik loglik;
	
	struct bfgs opt;
	enum bfgs_task task;
	double f;
	struct vector grad;
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


#endif /* _IPROC_FIT_H */
