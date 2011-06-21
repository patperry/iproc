#include "port.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include "linesearch.h"

#define TASK_START	"START"
#define TASK_FG		"FG"
#define TASK_CONV	"CONV"
#define TASK_WARN	"WARN"
#define TASK_ERROR	"ERROR"

extern int dcsrch_(double *stp, double *f, double *g,
		   double *ftol, double *gtol, double *xtol, char *task,
		   double *stpmin, double *stpmax, f77int *isave,
		   double *dsave, ssize_t task_len);

static void linesearch_dcscrch(struct linesearch *ls)
{
	assert(ls);
	assert(linesearch_ctrl_valid(&ls->ctrl));

	ssize_t task_len = sizeof(ls->task) - 1;
	dcsrch_(&ls->stp, &ls->f, &ls->g, &ls->ctrl.ftol,
		&ls->ctrl.gtol, &ls->ctrl.xtol,
		ls->task, &ls->ctrl.stpmin, &ls->ctrl.stpmax,
		ls->isave, ls->dsave, task_len);
	assert(task_len == 60);

	assert(strncmp(ls->task, TASK_ERROR, strlen(TASK_ERROR)) != 0);
}

void linesearch_start(struct linesearch *ls, double f0, double g0,
		      const struct linesearch_ctrl *ctrl)
{
	assert(ls);
	assert(!isnan(f0));
	assert(g0 < 0.0);
	assert(ctrl);
	assert(linesearch_ctrl_valid(ctrl));

	ls->ctrl = *ctrl;
	ls->stp = NAN;
	ls->f = f0;
	ls->g = g0;
	ls->done = false;
}

double linesearch_advance(struct linesearch *ls, double stp, double f, double g)
{
	assert(ls);
	assert(linesearch_ctrl_valid(&ls->ctrl));
	assert(stp > 0);
	assert(isfinite(stp));
	assert(isfinite(f));
	assert(isfinite(g));
	assert(ls->ctrl.stpmin <= stp && stp <= ls->ctrl.stpmax);
	assert(isnan(ls->stp) || stp == ls->stp);
	assert(isnan(ls->stp)
	       || strncmp(ls->task, TASK_FG, strlen(TASK_FG)) == 0);

	if (isnan(ls->stp)) {
		ls->stp = stp;

		strcpy(ls->task, TASK_START);
		linesearch_dcscrch(ls);
		linesearch_advance(ls, stp, f, g);
	}

	ls->f = f;
	ls->g = g;
	linesearch_dcscrch(ls);

	if (strncmp(ls->task, TASK_FG, strlen(TASK_FG)) == 0) {
		return ls->stp;
	} else {
		assert(strncmp(ls->task, TASK_CONV, strlen(TASK_CONV)) == 0
		       || strncmp(ls->task, TASK_WARN, strlen(TASK_WARN)) == 0);
		ls->done = true;
		return 0;
	}

	return ls->stp;
}

bool linesearch_converged(const struct linesearch *ls)
{
	assert(ls);
	if (!ls->done) {
		return false;
	} else {
		return (strncmp(ls->task, TASK_CONV, strlen(TASK_CONV)) == 0);
	}
}

const char *linesearch_error(const struct linesearch *ls)
{
	assert(ls);
	assert(ls->done);

	if (!ls->done) {
		return "LINESEARCH FAILED TO CONVERGE";
	} else {
		return ls->task + strlen("WARNING: ");
	}
}
