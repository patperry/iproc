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

struct recv_fit_constr {
	struct matrix ce;
	struct vector be;
};

struct recv_fit_resid {
	struct vector vector;
	double norm2;
};

struct recv_fit_eval {
	struct vector params;
	struct matrix coefs;
	struct vector duals;
	struct recv_loglik loglik;
	struct recv_fit_resid resid;
	bool in_domain;
};

struct recv_fit_kkt {
	struct matrix matrix;
	struct ldlfac ldl;
	enum matrix_uplo uplo;
	bool factored;
};

struct recv_fit_search {
	struct vector vector;
};

struct recv_fit_rgrad {
	struct vector vector;
};

struct recv_fit {
	struct recv_fit_ctrl ctrl;
	const struct messages *msgs;
	const struct design *design;
	const struct actors *senders;

	/* working frame + model */
	struct frame frame;
	struct recv_model model;

	/* optimization constraints */
	struct recv_fit_constr constr;

	/* optimization workspace */
	struct recv_fit_eval eval[2], *cur, *prev;
	struct recv_fit_kkt kkt;
	struct recv_fit_search search;
	struct recv_fit_rgrad rgrad;

	/* additional workspace */
	struct linesearch ls;
	enum recv_fit_task task;
	double step;
};

void recv_fit_init(struct recv_fit *fit,
		   const struct messages *msgs,
		   const struct design *design,
		   const struct actors *senders,
		   const struct matrix *coefs0,
		   const struct recv_fit_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);

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
