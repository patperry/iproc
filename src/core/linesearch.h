#ifndef _LINESEARCH_H
#define _LINESEARCH_H

#include <float.h>

#define LINSEARCH_FTOL0		1e-4
#define LINSEARCH_GTOL0		0.9
#define LINSEARCH_XTOL0		DBL_EPSILON
#define LINSEARCH_STPMIN0	DBL_EPSILON
#define LINSEARCH_STPMAX0	(1.0/DBL_EPSILON)


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

/* default control */
const struct linesearch_ctrl LINSEARCH_CTRL0 = {
	LINSEARCH_FTOL0,
	LINSEARCH_GTOL0,
	LINSEARCH_XTOL0,
	LINSEARCH_STPMIN0,
	LINSEARCH_STPMAX0
};

/* start/advance */
void linesearch_start(struct linesearch *ls, double f0, double g0,
		      const struct linesearch_ctrl *ctrl);
double linesearch_advance(struct linesearch *ls, double stp, double f, double g);

/* convergence failure */
const char *linesearch_warning(const struct linesearch *ls);

/* optimal values */
static inline double linesearch_step(const struct linesearch *ls);
static inline double linesearch_value(const struct linesearch *ls);
static inline double linesearch_grad(const struct linesearch *ls);

/* control parameters */
static inline bool linesearch_control_valid(const struct linesearch_ctrl *ctrl);


/* inline function definitions */
bool linesearch_control_valid(const struct linesearch_ctrl *ctrl)
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
	assert(ls->done);
	return ls->stp;
}

double linesearch_value(const struct linesearch *ls)
{
	assert(ls);
	assert(ls->done);
	return ls->f;
}

double linesearch_grad(const struct linesearch *ls)
{
	assert(ls);
	assert(ls->done);
	return ls->g;
}

#endif /* _LINESEARCH_H */
