#include "port.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include "linesearch.h"

#define DEFAULT_FTOL 1e-4
#define DEFAULT_GTOL 0.9
#define DEFAULT_XTOL 0.1

#define TASK_START	"START"
#define TASK_FG		"FG"
#define TASK_CONV	"CONV"
#define TASK_WARN	"WARN"
#define TASK_ERROR	"ERROR"

extern int dcsrch_(double *stp, double *f, double *g,
		   double *ftol, double *gtol, double *xtol, char *task,
		   double *stpmin, double *stpmax, f77int *isave,
		   double *dsave, f77int task_len);


static void linesearch_dcscrch(struct linesearch *ls)
{
	f77int task_len = sizeof(ls->task) - 1;
	dcsrch_(&ls->stp, &ls->f, &ls->g, &ls->ftol, &ls->gtol, &ls->xtol,
		ls->task, &ls->stpmin, &ls->stpmax, ls->isave, ls->dsave,
		task_len);
	assert(strncmp(ls->task, TASK_ERROR, strlen(TASK_ERROR)) != 0);
}

static void linesearch_init(struct linesearch *ls, double f0, double g0)
{
	assert(ls);
	assert(!isnan(f0));
	assert(g0 < 0.0);
	
	ls->stp = NAN;	
	ls->f = f0;
	ls->g = g0;
	ls->done = false;
	linesearch_set_fgtol(ls, DEFAULT_FTOL, DEFAULT_GTOL);
	linesearch_set_xtol(ls, DEFAULT_XTOL);
	linesearch_set_stpminmax(ls, 0, -f0 / (DEFAULT_GTOL * g0));
}

struct linesearch linesearch_make(double f0, double g0)
{
	struct linesearch ls;
	linesearch_init(&ls, f0, g0);
	return ls;
}

double linesearch_advance(struct linesearch *ls, double stp, double f, double g)
{
	assert(ls);
	assert(stp > 0);
	assert(isfinite(stp));
	assert(isfinite(f));
	assert(isfinite(g));
	assert(linesearch_stpmin(ls) <= stp && stp <= linesearch_stpmax(ls));
	assert(isnan(ls->stp) || stp == ls->stp);
	assert(isnan(ls->stp) || strncmp(ls->task, TASK_FG, strlen(TASK_FG)) == 0);	

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

const char *linesearch_warning(const struct linesearch *ls)
{
	assert(ls);
	assert(ls->done);
	
	if (strncmp(ls->task, TASK_CONV, strlen(TASK_CONV)) == 0) {
		return NULL;
	} else {
		return ls->task + strlen("WARNING: ");
	}
}

