#ifndef _LINESEARCH_H
#define _LINESEARCH_H

struct linesearch {
	double stp, f, g, ftol, gtol, xtol, stpmin, stpmax;
	char task[61];
	f77int isave[2];
	double dsave[13];
	bool done;
};

/* init/advance */
struct linesearch linesearch_make(double f0, double g0);
double linesearch_advance(struct linesearch *ls, double stp, double f, double g);

/* convergence failure */
const char *linesearch_warning(const struct linesearch *ls);

/* optimal values */
static inline double linesearch_step(const struct linesearch *ls);
static inline double linesearch_value(const struct linesearch *ls);
static inline double linesearch_grad(const struct linesearch *ls);

/* control parameters */
static inline double linesearch_ftol(const struct linesearch *ls);
static inline double linesearch_gtol(const struct linesearch *ls);
static inline void linesearch_set_fgtol(struct linesearch *ls, double ftol, double gtol);
static inline double linesearch_xtol(const struct linesearch *ls);
static inline void linesearch_set_xtol(struct linesearch *ls, double xtol);
static inline double linesearch_stpmin(const struct linesearch *ls);
static inline double linesearch_stpmax(const struct linesearch *ls);
static inline void linesearch_set_stpminmax(struct linesearch *ls, double stpmin, double stpmax);


/* inline function definitions */
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

double linesearch_ftol(const struct linesearch *ls)
{
	assert(ls);
	return ls->ftol;
}

double linesearch_gtol(const struct linesearch *ls)
{
	assert(ls);
	return ls->gtol;
}

void linesearch_set_fgtol(struct linesearch *ls, double ftol, double gtol)
{
	assert(ls);
	assert(ftol >= 0.0);
	assert(gtol >= 0.0);
	assert(ftol < gtol);
	ls->ftol = ftol;
	ls->gtol = gtol;
}

double linesearch_xtol(const struct linesearch *ls)
{
	assert(ls);
	return ls->xtol;
}

void linesearch_set_xtol(struct linesearch *ls, double xtol)
{
	assert(ls);
	assert(xtol >= 0.0); 
	ls->xtol = xtol;
}

double linesearch_stpmin(const struct linesearch *ls)
{
	assert(ls);
	return ls->stpmin;
}

double linesearch_stpmax(const struct linesearch *ls)
{
	assert(ls);
	return ls->stpmax;
}

void linesearch_set_stpminmax(struct linesearch *ls, double stpmin, double stpmax)
{
	assert(ls);
	assert(stpmin >= 0.0);
	assert(stpmax >= stpmin);
	ls->stpmin = stpmin;
	ls->stpmax = stpmax;
}

#endif /* _LINESEARCH_H */
