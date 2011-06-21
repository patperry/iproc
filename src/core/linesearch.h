#ifndef _LINESEARCH_H
#define _LINESEARCH_H

#include <float.h>

#define LINESEARCH_FTOL0	(1e-4)
#define LINESEARCH_GTOL0	(0.9)
#define LINESEARCH_XTOL0	(1e-15)
#define LINESEARCH_STPMIN0	(1e-15)
#define LINESEARCH_STPMAX0	(1e+15)

#define LINESEARCH_CTRL0 { \
		LINESEARCH_FTOL0, \
		LINESEARCH_GTOL0, \
		LINESEARCH_XTOL0, \
		LINESEARCH_STPMIN0, \
		LINESEARCH_STPMAX0 \
	}


struct linesearch_ctrl {
	double ftol;
	double gtol;
	double xtol;
	double stpmin;
	double stpmax;
};

struct linesearch {
	struct linesearch_ctrl ctrl;
	double stp, f, g;
	char task[61];
	f77int isave[2];
	double dsave[13];
	bool done;
};

/* start/advance */
void linesearch_start(struct linesearch *ls, double f0, double g0,
		      const struct linesearch_ctrl *ctrl);
double linesearch_advance(struct linesearch *ls, double stp, double f, double g);

/* convergence/error */
bool linesearch_converged(const struct linesearch *ls);
const char *linesearch_error(const struct linesearch *ls);

/* optimal values */
static inline double linesearch_step(const struct linesearch *ls);
static inline double linesearch_value(const struct linesearch *ls);
static inline double linesearch_grad(const struct linesearch *ls);

/* control parameters */
static inline bool linesearch_ctrl_valid(const struct linesearch_ctrl *ctrl);


/* inline function definitions */
bool linesearch_ctrl_valid(const struct linesearch_ctrl *ctrl)
{
	assert(ctrl);

	if (!(ctrl->ftol >= 0)) {
		return false;
	} else if (!(ctrl->gtol >= 0)) {
		return false;
	} else if (!(ctrl->xtol >= 0)) {
		return false;
	} else if (!(ctrl->stpmin >= 0)) {
		return false;
	} else if (!(ctrl->stpmax >= ctrl->stpmin)) {
		return false;
	} else {
		return true;
	}
}

double linesearch_step(const struct linesearch *ls)
{
	assert(ls);
	return ls->stp;
}

double linesearch_value(const struct linesearch *ls)
{
	assert(ls);
	return ls->f;
}

double linesearch_grad(const struct linesearch *ls)
{
	assert(ls);
	return ls->g;
}

#endif /* _LINESEARCH_H */
