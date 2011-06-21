#ifndef _BFGS_H
#define _BFGS_H

#include "linesearch.h"
#include "matrix.h"
#include "vector.h"

#define BFGS_ABSTOL0	(0.1)
#define BFGS_RELTOL0	(1e-5)
#define BFGS_LSMAX0	(20)
#define BFGS_LSCTRL0	LINESEARCH_CTRL0

#define BFGS_CTRL0 { \
		BFGS_ABSTOL0, \
		BFGS_RELTOL0, \
		BFGS_LSMAX0, \
		BFGS_LSCTRL0 \
	}


struct bfgs_ctrl {
	double abstol;	
	double reltol;
	ssize_t ls_maxit;	
	struct linesearch_ctrl ls;
};

struct bfgs {
	struct bfgs_ctrl ctrl;
	struct matrix inv_hess;
	struct vector search_dir;
	struct vector grad0;
	struct vector x0;
	struct vector x;
	double f0;
	double g0;
	double step;
	ssize_t ls_it;
	bool done;
	const char *errmsg;

	/* workspace */
	struct linesearch ls;
	struct vector s;
	struct vector y;
	struct vector H_y;
	bool first_step;
};



void bfgs_init(struct bfgs *opt, ssize_t n, const struct bfgs_ctrl *ctrl);
void bfgs_deinit(struct bfgs *opt);

static inline ssize_t bfgs_dim(const struct bfgs *opt);

const struct vector *bfgs_start(struct bfgs *opt, const struct vector *x0,
				double f0, const struct vector *grad0);
const struct vector *bfgs_advance(struct bfgs *opt, double f, const struct vector *grad);


/* convergence/warning */
bool bfgs_converged(const struct bfgs *opt);
const char *bfgs_error(const struct bfgs *opt);

/* current values */
static inline const struct vector *bfgs_location(const struct bfgs *opt);
static inline double bfgs_value(const struct bfgs *opt);
static inline const struct vector *bfgs_grad(const struct bfgs *opt);
static inline const struct matrix *bfgs_inv_hess(const struct bfgs *opt);
static inline double bfgs_decrement(const struct bfgs *opt);

/* control parameters */
static inline bool bfgs_ctrl_valid(const struct bfgs_ctrl *ctrl);


/* inline function definitions */
bool bfgs_ctrl_valid(const struct bfgs_ctrl *ctrl)
{
	assert(ctrl);
	
	return linesearch_ctrl_valid(&ctrl->ls);
}

ssize_t bfgs_dim(const struct bfgs *opt)
{
	assert(opt);
	return vector_dim(&opt->x0);
}

const struct vector *bfgs_location(const struct bfgs *opt)
{
	assert(opt);
	assert(!opt->errmsg);
	return &opt->x0;
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
	assert(!opt->first_step);
	return &opt->inv_hess;
}

double bfgs_decrement(const struct bfgs *opt)
{
	assert(opt);
	return opt->g0;
}

#endif /* _BFGS_H */
