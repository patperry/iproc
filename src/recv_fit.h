#ifndef RECV_FIT_H
#define RECV_FIT_H

#include "constr.h"
#include "design.h"
#include "linesearch.h"
#include "messages.h"

#include "recv_loglik.h"
#include "recv_model.h"

#define RECV_FIT_GTOL0		(1e-8)
#define RECV_FIT_XTOL0		(1e-8)
#define RECV_FIT_LSMAX0		(10)
#define RECV_FIT_LSCTRL0	LINESEARCH_CTRL0

#define RECV_FIT_CTRL0 ((struct recv_fit_ctrl) { \
		RECV_FIT_GTOL0, \
		RECV_FIT_XTOL0, \
		RECV_FIT_LSMAX0, \
		RECV_FIT_LSCTRL0 \
	})

struct recv_fit_ctrl {
	double gtol;
	double xtol;
	size_t ls_maxit;
	struct linesearch_ctrl ls;
};

enum recv_fit_task {
	RECV_FIT_CONV = 0,
	RECV_FIT_STEP = 1,
	RECV_FIT_ERR_LNSRCH = -1,	// linesearch failed to converge
	RECV_FIT_ERR_XTOL = -2,	// step size is smaller than tolerance
};


struct recv_params {
	double *all;
	size_t dim;
	struct recv_coefs coefs;
	double *duals;
	int owner;
};

struct recv_fit_resid {
	struct recv_params params;
	double norm2;
};

struct recv_fit_eval {
	struct recv_loglik loglik;
	struct recv_params params;
	struct recv_fit_resid resid;
};

struct recv_fit_kkt {
	double *matrix;
	ptrdiff_t *ldl_ipiv;
	double *ldl_work;
	double *imat_buf;
	size_t ldl_lwork;
	int factored;
};

struct recv_fit_search {
	struct recv_params params;
};

struct recv_fit_rgrad {
	struct recv_params params;
};

struct recv_fit {
	struct recv_fit_ctrl ctrl;
	const struct messages *xmsgs, *ymsgs;

	/* working frame + model */
	struct frame *frame;
	struct recv_model model;

	/* null deviance */
	double dev0;

	/* optimization constraints */
	struct constr *constr;

	/* optimization workspace */
	struct recv_fit_eval eval[2], *cur, *prev;
	struct recv_fit_kkt kkt;
	struct recv_fit_search search;
	struct recv_fit_rgrad rgrad;

	/* additional workspace */
	//struct linesearch ls;
	enum recv_fit_task task;
	double step;
};

void recv_fit_init(struct recv_fit *fit,
		   struct frame *f,
		   struct constr *c,
		   const struct messages *xmsgs,
		   const struct messages *ymsgs,
		   const struct recv_fit_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);


void recv_params_init(struct recv_params *params, const struct frame *f,
		      const struct constr *c);
void recv_params_init_view(struct recv_params *params, const struct frame *f,
			   const struct constr *c, const double *data);
void recv_params_deinit(struct recv_params *params);


/* fitting */
enum recv_fit_task recv_fit_start(struct recv_fit *fit, const struct recv_params *params0);
enum recv_fit_task recv_fit_advance(struct recv_fit *fit);
const char *recv_fit_errmsg(const struct recv_fit *fit);

/* current values */
const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit);
double recv_fit_dev(const struct recv_fit *fit);
double recv_fit_dev0(const struct recv_fit *fit);
const struct recv_coefs *recv_fit_coefs(const struct recv_fit *fit);
const double *recv_fit_duals(const struct recv_fit *fit);
double recv_fit_step(const struct recv_fit *fit);
double recv_fit_grad_norm2(const struct recv_fit *fit);

/* control parameters */
static inline int recv_fit_ctrl_valid(const struct recv_fit_ctrl *ctrl);

/* inline function definitions */
int recv_fit_ctrl_valid(const struct recv_fit_ctrl *ctrl)
{
	assert(ctrl);

	if (!(ctrl->gtol > 0)) {
		return 0;
	} else if (!(ctrl->ls_maxit > 0)) {
		return 0;
	} else {
		return linesearch_ctrl_valid(&ctrl->ls);
	}
}

#endif /* RECV_FIT_H */
