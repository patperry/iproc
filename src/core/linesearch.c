#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include "f77.h"

#include "linesearch.h"

#define TASK_START		"START"
#define TASK_FG			"FG"
#define TASK_CONV		"CONV"
#define TASK_WARN		"WARN"
#define TASK_WARN_ROUND		"WARNING: ROUNDING ERRORS PREVENT PROGRESS"
#define TASK_WARN_XTOL		"WARNING: XTOL TEST SATISFIED"
#define TASK_WARN_STPMAX	"WARNING: STP = STPMAX"
#define TASK_WARN_STPMIN	"WARNING: STP = STPMIN"
#define TASK_ERROR		"ERROR"

extern int dcsrch_(double *stp, double *f, double *g,
		   double *ftol, double *gtol, double *xtol, char *task,
		   double *stpmin, double *stpmax, f77int *isave,
		   double *dsave, size_t task_len);

static void linesearch_dcsrch(struct linesearch *ls)
{
	assert(ls);
	assert(isfinite(ls->f));
	assert(isfinite(ls->g));
	assert(linesearch_ctrl_valid(&ls->ctrl));
	assert((ls->task == LINESEARCH_STEP
		&& !strncmp(ls->taskbuf, TASK_FG, strlen(TASK_FG)))
	       || !strncmp(ls->taskbuf, TASK_START, strlen(TASK_START)));

	//if (ls->task != LINESEARCH_STEP) {
	//      fprintf(stderr, "LNSRCH: f0 = %.22e, g0 = %.22e\n", ls->f, ls->g);
	//} else {
	//      fprintf(stderr, "LNSRCH: stp = %.22e, f = %.22e, g = %.22e\n", ls->stp, ls->f, ls->g);
	//}

	size_t task_len = sizeof(ls->taskbuf) - 1;

	assert(task_len == 60);
	dcsrch_(&ls->stp, &ls->f, &ls->g, &ls->ctrl.ftol,
		&ls->ctrl.gtol, &ls->ctrl.xtol,
		ls->taskbuf, &ls->ctrl.stpmin, &ls->ctrl.stpmax,
		(f77int *)ls->isave, ls->dsave, task_len);

	if (!strncmp(ls->taskbuf, TASK_CONV, strlen(TASK_CONV))) {
		ls->task = LINESEARCH_CONV;
	} else if (!strncmp(ls->taskbuf, TASK_FG, strlen(TASK_FG))) {
		ls->task = LINESEARCH_STEP;
	} else
	    if (!strncmp(ls->taskbuf, TASK_WARN_ROUND, strlen(TASK_WARN_ROUND)))
	{
		ls->task = LINESEARCH_WARN_ROUND;
	} else
	    if (!strncmp(ls->taskbuf, TASK_WARN_XTOL, strlen(TASK_WARN_XTOL))) {
		ls->task = LINESEARCH_WARN_XTOL;
	} else
	    if (!strncmp
		(ls->taskbuf, TASK_WARN_STPMAX, strlen(TASK_WARN_STPMAX))) {
		ls->task = LINESEARCH_WARN_STPMAX;
	} else
	    if (!strncmp
		(ls->taskbuf, TASK_WARN_STPMIN, strlen(TASK_WARN_STPMIN))) {
		ls->task = LINESEARCH_WARN_STPMIN;
	} else {
		assert(!strncmp(ls->taskbuf, TASK_ERROR, sizeof(TASK_ERROR)));
		assert(0);
	}

	//if (ls->task < 0) {
	//      fprintf(stderr, "[LINESEARCH] %s\n", linesearch_warnmsg(ls->task));
	//}
}

void linesearch_start(struct linesearch *ls, double stp0, double f0, double g0,
		      const struct linesearch_ctrl *ctrl)
{
	assert(ls);
	assert(isfinite(stp0));
	assert(isfinite(f0));
	assert(g0 < 0.0);
	assert(ctrl);
	assert(linesearch_ctrl_valid(ctrl));
	assert(ctrl->stpmin <= stp0 && stp0 <= ctrl->stpmax);

	ls->ctrl = *ctrl;
	ls->stp = stp0;
	ls->f = f0;
	ls->f0 = f0;
	ls->g = g0;
	ls->g0 = g0;

	strcpy(ls->taskbuf, TASK_START);
	linesearch_dcsrch(ls);
}

enum linesearch_task linesearch_advance(struct linesearch *ls, double f,
					double g)
{
	assert(ls);
	assert(linesearch_ctrl_valid(&ls->ctrl));
	assert(isfinite(f));
	assert(isfinite(g));
	assert(ls->task = LINESEARCH_STEP);

	ls->f = f;
	ls->g = g;
	linesearch_dcsrch(ls);
	return ls->task;
}

const char *linesearch_warnmsg(enum linesearch_task task)
{
	switch (task) {
	case LINESEARCH_STEP:
		return "search in progress";
	case LINESEARCH_WARN_ROUND:
		return "rounding errors prevent progress";
	case LINESEARCH_WARN_XTOL:
		return "xtol test satisfied";
	case LINESEARCH_WARN_STPMAX:
		return "stp = stpmax";
	case LINESEARCH_WARN_STPMIN:
		return "stp = stpmin";
	case LINESEARCH_CONV:
		return NULL;
	}
	assert(0);
	return NULL;
}
