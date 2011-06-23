#ifndef _BFGS_H
#define _BFGS_H

#include "linesearch.h"
#include "matrix.h"
#include "vector.h"

#define BFGS_GTOL0	(1e-5)
#define BFGS_LSMAX0	(20)
#define BFGS_LSCTRL0	LINESEARCH_CTRL0

#define BFGS_CTRL0 { \
		BFGS_GTOL0, \
		BFGS_LSMAX0, \
		BFGS_LSCTRL0 \
	}

struct bfgs_ctrl {
	double gtol;
	ssize_t ls_maxit;
	struct linesearch_ctrl ls;
};

enum bfgs_task {
	BFGS_CONV = 0,
	BFGS_STEP = 1,
	BFGS_ERR_LNSRCH = -1, // linesearch failed to converge
	BFGS_OVFLW_GRAD = -2, // overflow in computing norm of gradient
};

struct bfgs {
	/* control/status */
	struct bfgs_ctrl ctrl;
	enum bfgs_task task;
	bool first_step;

	/* current values */
	double f0;
	struct vector grad0;
	struct vector x0;
	
	/* next step */
	struct vector x;
	struct vector search;	
	struct vector step;
	
	/* linsearch workspace */
	struct linesearch ls;
	ssize_t ls_it;

	/* inverse hessian estimate */
	struct matrix inv_hess;

	/* auxiliary variables */
	double df;
	struct vector dg;
	struct vector H_dg;
};

void bfgs_init(struct bfgs *opt, ssize_t n, const struct bfgs_ctrl *ctrl);
void bfgs_deinit(struct bfgs *opt);

static inline ssize_t bfgs_dim(const struct bfgs *opt);

enum bfgs_task bfgs_start(struct bfgs *opt, const struct vector *x0,
			  double f0, const struct vector *grad0);
static inline const struct vector *bfgs_next(const struct bfgs *opt);
enum bfgs_task bfgs_advance(struct bfgs *opt, double f,
			    const struct vector *grad);
const char *bfgs_errmsg(enum bfgs_task task);


/* current values */
static inline const struct vector *bfgs_cur(const struct bfgs *opt);
static inline double bfgs_value(const struct bfgs *opt);
static inline const struct vector *bfgs_grad(const struct bfgs *opt);
static inline const struct matrix *bfgs_inv_hess(const struct bfgs *opt);

/* control parameters */
static inline bool bfgs_ctrl_valid(const struct bfgs_ctrl *ctrl);

/* inline function definitions */
bool bfgs_ctrl_valid(const struct bfgs_ctrl *ctrl)
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

ssize_t bfgs_dim(const struct bfgs *opt)
{
	assert(opt);
	return vector_dim(&opt->x0);
}

const struct vector *bfgs_cur(const struct bfgs *opt)
{
	assert(opt);
	return &opt->x0;
}

const struct vector *bfgs_next(const struct bfgs *opt)
{
	assert(opt);
	assert(opt->task == BFGS_STEP);
	return &opt->x;
}

double bfgs_value(const struct bfgs *opt)
{
	assert(opt);
	return opt->f0;
}

const struct vector *bfgs_grad(const struct bfgs *opt)
{
	assert(opt);
	return &opt->grad0;
}

const struct matrix *bfgs_inv_hess(const struct bfgs *opt)
{
	assert(opt);
	if (opt->first_step) {
		return NULL;
	} else {
		return &opt->inv_hess;
	}
}

#endif /* _BFGS_H */
