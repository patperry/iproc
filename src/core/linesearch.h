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

/*
p <= stmin || *stp >= stmax)) {
	s_copy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS", task_len, (
									     ftnlen)41);
}
if (brackt && stmax - stmin <= *xtol * stmax) {
	s_copy(task, "WARNING: XTOL TEST SATISFIED", task_len, (ftnlen)28);
}
if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
	s_copy(task, "WARNING: STP = STPMAX", task_len, (ftnlen)21);
}
if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
	s_copy(task, "WARNING: STP = STPMIN", task_len, (ftnlen)21);
}
if (*f <= ftest && abs(*g) <= *gtol * (-ginit)) {
	s_copy(task, "CONVERGENCE
*/
	       
enum linesearch_task {
	LINESEARCH_CONV = 0,
	LINESEARCH_STEP = 1,
	LINESEARCH_WARN_ROUND = -1,
	LINESEARCH_WARN_XTOL = -2,
	LINESEARCH_WARN_STPMAX = -3,
	LINESEARCH_WARN_STPMIN = -4,
};

struct linesearch {
	struct linesearch_ctrl ctrl;
	double stp, f, f0, g, g0;
	enum linesearch_task task;
	char taskbuf[61];
	f77int isave[2];
	double dsave[13];
};

/* start/advance */
void linesearch_start(struct linesearch *ls, double stp0, double f0, double g0,
		      const struct linesearch_ctrl *ctrl);
enum linesearch_task linesearch_advance(struct linesearch *ls, double f, double g);

/* weak Wolfe conditions */
static inline bool linesearch_sdec(const struct linesearch *ls);  // sufficient decrease (Armijo)
static inline bool linesearch_curv(const struct linesearch *ls);  // curvature

/* convergence/error */
const char *linesearch_warnmsg(enum linesearch_task task);


/* optimal values */
static inline double linesearch_step(const struct linesearch *ls);
static inline double linesearch_value(const struct linesearch *ls);
static inline double linesearch_grad(const struct linesearch *ls);

/* initial values */
static inline double linesearch_value0(const struct linesearch *ls);
static inline double linesearch_grad0(const struct linesearch *ls);

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

double linesearch_value0(const struct linesearch *ls)
{
	assert(ls);
	return ls->f0;
}

double linesearch_grad0(const struct linesearch *ls)
{
	assert(ls);
	return ls->g0;
}

double linesearch_value(const struct linesearch *ls)
{
	assert(ls);
	return ls->f;
}

bool linesearch_sdec(const struct linesearch *ls)
{
	assert(ls);
	double stp = linesearch_step(ls);
	double f0 = linesearch_value0(ls);
	double f = linesearch_value(ls);
	double g0 = linesearch_grad0(ls);
	double c1 = ls->ctrl.ftol;
	
	return f <= f0 + c1 * stp * g0;
}

bool linesearch_curv(const struct linesearch *ls)
{
	assert(ls);
	double g0 = linesearch_grad0(ls);
	double g = linesearch_grad(ls);
	double c2 = ls->ctrl.ftol;

	return g >= c2 * g0;
}

double linesearch_grad(const struct linesearch *ls)
{
	assert(ls);
	return ls->g;
}

#endif /* _LINESEARCH_H */
