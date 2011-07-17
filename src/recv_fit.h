#ifndef _RECV_FIT_H
#define _RECV_FIT_H

#include "bfgs.h"
#include "design.h"
#include "linalg.h"
#include "linesearch.h"
#include "matrix.h"
#include "messages.h"
#include "recv_loglik.h"
#include "recv_model.h"
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
	RECV_FIT_ERR_LNSRCH = -1,	// linesearch failed to converge
	RECV_FIT_ERR_XTOL = -2,	// step size is smaller than tolerance
};


struct recv_fit_cohort {
	/* null log-likelihood */
	struct matrix imat0;
	struct vector score0;
	double dev0;
	
	/* regularization terms */
	struct vector scale;	/* sample variances of the covariates */

	/* optimization problem */
	struct matrix ce;	/* equality constraints: ce' * coef = be */
	struct vector be;	/* cont'd */
	ssize_t ne;
	
	/* optimization problem workspace */
	struct vector params0;	/* initial primal and dual parameters */	
	struct vector params;	/* current primal and dual parameters */
	struct vector resid;	/* current dual and primal residuals */
	struct matrix kkt;	/* current KKT matrix */

	/* linesearch */
	struct vector search;   /* search direction */
	double rss;		/* residual sum of squares */
	double grss;            /* directional derivative of rss */	
	struct vector grad_rss;	/* gradient of rss */
	
	/* additional workspace */
	struct ldlfac ldl;
};

struct recv_fit {
	struct recv_fit_ctrl ctrl;
	const struct design *design;
	const struct actors *senders;
	const struct messages *msgs;

	struct frame frame;
	struct matrix coefs;
	struct recv_model model;
	struct recv_loglik loglik;
	struct recv_fit_cohort *cohorts;
		
	/* regularization */
	double penalty;
	
	/* additional workspace */
	struct linesearch ls;
	struct symeig eig;

	enum recv_fit_task task;
	double step;
	double rss;
	double grss;
};

void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct actors *senders,
		   const struct matrix *coefs0,
		   const struct recv_fit_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);

/* problem constraints */
ssize_t recv_fit_rank(const struct recv_fit *fit, ssize_t c);
const struct matrix *recv_fit_ce(const struct recv_fit *fit, ssize_t c);
const struct vector *recv_fit_be(const struct recv_fit *fit, ssize_t c);

/* null values */
double recv_fit_dev0(const struct recv_fit *fit, ssize_t c);
const struct vector *recv_fit_score0(const struct recv_fit *fit, ssize_t c);
const struct matrix *recv_fit_imat0(const struct recv_fit *fit, ssize_t c);

enum recv_fit_task recv_fit_advance(struct recv_fit *fit);
const char *recv_fit_errmsg(const struct recv_fit *fit);

/* current values */
const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit);
const struct matrix *recv_fit_coefs(const struct recv_fit *fit);
double recv_fit_step(const struct recv_fit *fit);
double recv_fit_grad_norm2(const struct recv_fit *fit);

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
