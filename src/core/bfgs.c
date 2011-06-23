#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "bfgs.h"

static void linesearch(struct bfgs *opt, double step0, double f0, double g0)
{
	opt->ls_it = 0;
	vector_assign_copy(&opt->x, &opt->x0);
	vector_axpy(step0, &opt->search, &opt->x);

	linesearch_start(&opt->ls, step0, f0, g0, &opt->ctrl.ls);
}

static void update(struct bfgs *opt, double f, const struct vector *grad)
{
	struct matrix *H = &opt->inv_hess;
	struct vector *s = &opt->step;
	struct vector *y = &opt->dg;

	/* update step */
	vector_assign_copy(s, &opt->search);
	vector_scale(s, linesearch_step(&opt->ls));

	/* update df */
	opt->df = f - opt->f0;

	/* update dg */
	vector_assign_copy(y, grad);
	vector_axpy(-1.0, &opt->grad0, y);

	double s_y = vector_dot(s, y);
	
	/* NOTE: could use damped update instead (Nocedal and Wright, p. 537) */	
	assert(s_y > 0); 

	/* initialize inv hessian on first step (Nocedal and Wright, p. 143) */
	if (opt->first_step) { /*  */
		double y_y = vector_dot(y, y);
		assert(y_y > 0);
		double scale = s_y / y_y;
		ssize_t i, n = vector_dim(y);

		matrix_fill(H, 0.0);
		for (i = 0; i < n; i++) {
			matrix_set_item(H, i, i, scale);
		}
		opt->first_step = false;
	}

	/* compute H_y */
	struct vector *H_y = &opt->H_dg;
	matrix_mul(1.0, TRANS_NOTRANS, H, y, 0.0, H_y);

	double y_H_y = vector_dot(H_y, y);
	double scale1 = (1.0 + (y_H_y / s_y)) / s_y;
	double rho = 1.0 / s_y;

	/* update inverse hessian */
	matrix_update1(H, scale1, s, s);
	matrix_update1(H, -rho, H_y, s);
	matrix_update1(H, -rho, s, H_y);
	
	/* update search direction */
	matrix_mul(-1.0, TRANS_NOTRANS, &opt->inv_hess, grad,
		   0.0, &opt->search);
	assert(isfinite(vector_norm(&opt->search)));

	/* update initial position, value, and grad */
	vector_assign_copy(&opt->x0, &opt->x);
	opt->f0 = f;
	vector_assign_copy(&opt->grad0, grad);
}

void bfgs_init(struct bfgs *opt, ssize_t n, const struct bfgs_ctrl *ctrl)
{
	assert(opt);
	assert(n >= 0);
	assert(ctrl);
	assert(bfgs_ctrl_valid(ctrl));

	opt->ctrl = *ctrl;
	matrix_init(&opt->inv_hess, n, n);
	vector_init(&opt->search, n);
	vector_init(&opt->grad0, n);
	vector_init(&opt->x0, n);
	vector_init(&opt->x, n);
	vector_init(&opt->step, n);
	vector_init(&opt->dg, n);
	vector_init(&opt->H_dg, n);
}

void bfgs_deinit(struct bfgs *opt)
{
	vector_deinit(&opt->H_dg);
	vector_deinit(&opt->dg);
	vector_deinit(&opt->step);
	vector_deinit(&opt->x);
	vector_deinit(&opt->x0);
	vector_deinit(&opt->grad0);
	vector_deinit(&opt->search);
	matrix_deinit(&opt->inv_hess);
}

enum bfgs_task bfgs_start(struct bfgs *opt, const struct vector *x0,
			  double f0, const struct vector *grad0)
{
	assert(opt);
	assert(x0);
	assert(isfinite(f0));
	assert(grad0);
	assert(vector_dim(x0) == bfgs_dim(opt));
	assert(vector_dim(grad0) == bfgs_dim(opt));

	double step0 = 1.0;

	opt->first_step = true;
	vector_assign_copy(&opt->x0, x0);
	opt->f0 = f0;
	vector_assign_copy(&opt->grad0, grad0);

#ifndef NDEBUG
	opt->df = NAN;
	vector_fill(&opt->step, NAN);
	vector_fill(&opt->dg, NAN);
	vector_fill(&opt->H_dg, NAN);	
	vector_fill(&opt->search, NAN);
	matrix_fill(&opt->inv_hess, NAN);
#endif
	
	double scale = vector_norm(grad0);

	assert(!isnan(scale));

	if (!isfinite(scale)) {
		opt->task = BFGS_OVFLW_GRAD;
	} else if (scale > 0) {
		opt->task = BFGS_STEP;

		vector_assign_copy(&opt->search, grad0);
		vector_scale(&opt->search, -1.0 / scale);
		double g0 = -scale;
		linesearch(opt, step0, f0, g0);
	} else {
		opt->task = BFGS_CONV;
	}
	return opt->task;
}

static bool converged(const struct bfgs *opt)
{
	double gtol = opt->ctrl.gtol;

	//fprintf(stderr, "|df| = %.22f; |grad| = %.22f; |step| = %.22f\n",
	//	fabs(opt->df /MAX(1, opt->f0)),
	//	vector_max_abs(&opt->grad0), vector_max_abs(&opt->step));
	
	ssize_t i, n = bfgs_dim(opt);
	for (i = 0; i < n; i++) {
		double x = vector_item(&opt->x0, i);
		double g = vector_item(&opt->grad0, i);
		
		if (!(fabs(g) < gtol * MAX(1.0, fabs(x))))
			return false;
	}
	return true;		
}

enum bfgs_task bfgs_advance(struct bfgs *opt, double f,
			    const struct vector *grad)
{
	assert(isfinite(f));
	assert(isfinite(vector_norm(grad)));
	assert(opt->task == BFGS_STEP);

	double g = vector_dot(grad, &opt->search);
	enum linesearch_task lstask = linesearch_advance(&opt->ls, f, g);
	bool ok = linesearch_sdec(&opt->ls) && linesearch_curv(&opt->ls);
	
	switch (lstask) {
	case LINESEARCH_CONV:
		break;

	case LINESEARCH_STEP:
		opt->ls_it++;

		if (opt->ls_it < opt->ctrl.ls_maxit) {
			vector_assign_copy(&opt->x, &opt->x0);
			vector_axpy(linesearch_step(&opt->ls), &opt->search, &opt->x);
			assert(opt->task == BFGS_STEP);
			goto out;
		} else if (ok) {
			break;
		} else {
			opt->task = BFGS_ERR_LNSRCH; // maximum number of iterations
		}
	default:
		if (ok) {
			break;
		} else {
			opt->task = BFGS_ERR_LNSRCH;
			goto out;
		}
	}

	update(opt, f, grad);

	// test for convergence
	if (converged(opt)) {
		opt->task = BFGS_CONV;
	} else {
		assert(opt->task == BFGS_STEP);
		
		double step0 = 1.0;
		double f0 = f;	
		double g0 = vector_dot(grad, &opt->search);
		assert(g0 < 0);
		
		linesearch(opt, step0, f0, g0);
	}
out:
	return opt->task;
}

const char *bfgs_errmsg(enum bfgs_task task)
{
	switch(task) {
		case BFGS_STEP:
			return "optimization in progress";
		case BFGS_ERR_LNSRCH:
			return "linesearch failed";
		case BFGS_OVFLW_GRAD:
			return "overflow computing norm of gradient";
		case BFGS_CONV:
			return NULL;
	}
	assert(0);
	return NULL;
}

