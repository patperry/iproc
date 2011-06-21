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

typedef struct _iproc_fit iproc_fit;

struct _iproc_fit {
	struct model *model;
	struct messages *messages;
	double penalty;
	struct recv_loglik *loglik;
	struct bfgs opt;
	double f;
	struct vector grad;
	const struct vector *xnext;
};

iproc_fit *iproc_fit_new(struct model *model0,
			 struct messages *messages, double penalty);
void iproc_fit_free(iproc_fit * fit);

bool iproc_fit_converged(iproc_fit * fit);
void iproc_fit_step(iproc_fit * fit);

#endif /* _IPROC_FIT_H */
