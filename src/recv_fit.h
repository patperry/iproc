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


#define RECV_FIT_GTOL0		(1e-8)
#define RECV_FIT_XTOL0		(1e-8)
#define RECV_FIT_LSMAX0		(10)
#define RECV_FIT_LSCTRL0	LINESEARCH_CTRL0
#define RECV_FIT_VARTOL		(1e-8)
#define RECV_FIT_EIGTOL		(1e-8)

#define RECV_FIT_CTRL0 ((struct recv_fit_ctrl) { \
		RECV_FIT_GTOL0, \
		RECV_FIT_XTOL0, \
		RECV_FIT_LSMAX0, \
		RECV_FIT_LSCTRL0, \
		RECV_FIT_VARTOL, \
		RECV_FIT_EIGTOL \
	})

struct recv_fit_ctrl {
	double gtol;
	double xtol;	
	ssize_t ls_maxit;
	struct linesearch_ctrl ls;
	double vartol;
	double eigtol;
};

enum recv_fit_task {
	RECV_FIT_CONV = 0,
	RECV_FIT_STEP = 1,
	RECV_FIT_ERR_LNSRCH = -1, // linesearch failed to converge
	RECV_FIT_ERR_XTOL = -2, // step size is smaller than tolerance
};


struct recv_fit {
	struct recv_fit_ctrl ctrl;
	const struct design *design;
	const struct messages *msgs;

	struct frame frame;
	struct vector coefs;
	struct model model;
	struct recv_loglik loglik;
	
	/* regularization terms */
	struct vector scale;  /* sample variances of the covariates */
	double penalty;	
	
	/* optimization problem */
	struct matrix ce_t;  /* equality constraints: ce * coef = be */
	struct vector be;    /* cont'd */
	ssize_t ne;

	/* optimization problem workspace */
	struct vector params; /* primal and dual parameters */
	struct vector resid;  /* dual and primal residuals */
	struct matrix kkt;    /* KKT matrix */
	struct vector search;
	double step;
	
	/* additional workspace */
	struct linesearch ls;
	struct symeig eig;	
	struct ldlfac ldl;
	enum recv_fit_task task;
};


void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct vector *coefs0,
		   const struct recv_fit_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);

/* problem constraints */

enum recv_fit_task recv_fit_advance(struct recv_fit *fit);
const char *recv_fit_errmsg(const struct recv_fit *fit);

/* current values */
const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit);
const struct vector *recv_fit_coefs(const struct recv_fit *fit);
double recv_fit_step(const struct recv_fit *fit);

/* control parameters */
static inline bool recv_fit_ctrl_valid(const struct recv_fit_ctrl *ctrl);

/* inline function definitions */
bool recv_fit_ctrl_valid(const struct recv_fit_ctrl *ctrl)
{
	assert(ctrl);
	
	if (!(ctrl->gtol > 0)) {
		return false;
	} else if (!(ctrl->ls_maxit > 0)) {
		return false;
	} else {
		return linesearch_ctrl_valid(&ctrl->ls);
	}
}


#endif /* _RECV_FIT_H */
