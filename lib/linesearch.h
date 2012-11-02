#ifndef LINESEARCH_H
#define LINESEARCH_H

#include <stddef.h>

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

enum linesearch_task {
	LINESEARCH_CONV = 0,
	LINESEARCH_STEP = 1,
	LINESEARCH_WARN_ROUND = -1,
	LINESEARCH_WARN_XTOL = -2,
	LINESEARCH_WARN_STPMAX = -3,
	LINESEARCH_WARN_STPMIN = -4,
};

/* control parameters */
static inline int linesearch_ctrl_valid(const struct linesearch_ctrl *ctrl);

/* inline function definitions */
int linesearch_ctrl_valid(const struct linesearch_ctrl *ctrl)
{
	assert(ctrl);

	if (!(ctrl->ftol >= 0)) {
		return 0;
	} else if (!(ctrl->gtol >= 0)) {
		return 0;
	} else if (!(ctrl->xtol >= 0)) {
		return 0;
	} else if (!(ctrl->stpmin >= 0)) {
		return 0;
	} else if (!(ctrl->stpmax >= ctrl->stpmin)) {
		return 0;
	} else {
		return 1;
	}
}

#endif /* LINESEARCH_H */
