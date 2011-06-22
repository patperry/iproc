#include "port.h"
#include <assert.h>
#include <math.h>
#include "bfgs.h"

static const struct vector *linesearch(struct bfgs *opt, double step0, double f0, double g0)
{
	opt->ls_it = 0;
	vector_assign_copy(&opt->x, &opt->x0);
	vector_axpy(step0, &opt->search_dir, &opt->x);

	linesearch_start(&opt->ls, step0, f0, g0, &opt->ctrl.ls);
	return &opt->x;
}

static void update_searchdir(struct bfgs *opt, const struct vector *grad)
{
	assert(opt);
	assert(!opt->first_step);

	matrix_mul(-1.0, TRANS_NOTRANS, &opt->inv_hess, grad,
		   0.0, &opt->search_dir);
	assert(isfinite(vector_norm(&opt->search_dir)));
}

static void update_hess(struct bfgs *opt, const struct vector *grad)
{
	struct matrix *H = &opt->inv_hess;
	struct vector *s = &opt->s;
	struct vector *y = &opt->y;

	/* compute s */
	vector_assign_copy(s, &opt->search_dir);
	vector_scale(s, linesearch_step(&opt->ls));

	/* compute y */
	vector_assign_copy(y, grad);
	vector_axpy(-1.0, &opt->grad0, y);

	double s_y = vector_dot(s, y);
	
	/* Note could use damped update instead (Nocedal and Wright, p. 537) */	
	assert(s_y > 0); 

	/* update H */
	if (opt->first_step) { /* Nocedal and Wright, p. 143 */
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
	struct vector *H_y = &opt->H_y;
	matrix_mul(1.0, TRANS_NOTRANS, H, y, 0.0, H_y);

	double y_H_y = vector_dot(H_y, y);
	double scale1 = (1.0 + (y_H_y / s_y)) / s_y;
	double rho = 1.0 / s_y;

	matrix_update1(H, scale1, s, s);
	matrix_update1(H, -rho, H_y, s);
	matrix_update1(H, -rho, s, H_y);
}

void bfgs_init(struct bfgs *opt, ssize_t n, const struct bfgs_ctrl *ctrl)
{
	assert(opt);
	assert(n >= 0);
	assert(ctrl);
	assert(bfgs_ctrl_valid(ctrl));

	opt->ctrl = *ctrl;
	matrix_init(&opt->inv_hess, n, n);
	vector_init(&opt->search_dir, n);
	vector_init(&opt->grad0, n);
	vector_init(&opt->x0, n);
	vector_init(&opt->x, n);
	vector_init(&opt->s, n);
	vector_init(&opt->y, n);
	vector_init(&opt->H_y, n);
}

void bfgs_deinit(struct bfgs *opt)
{
	vector_deinit(&opt->H_y);
	vector_deinit(&opt->y);
	vector_deinit(&opt->s);
	vector_deinit(&opt->x);
	vector_deinit(&opt->x0);
	vector_deinit(&opt->grad0);
	vector_deinit(&opt->search_dir);
	matrix_deinit(&opt->inv_hess);
}

const struct vector *bfgs_start(struct bfgs *opt, const struct vector *x0,
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
	vector_assign_copy(&opt->grad0, grad0);

	double scale = vector_norm(grad0);

	assert(!isnan(scale));

	if (!isfinite(scale)) {
		opt->done = true;
		opt->errmsg = "OVERFLOW IN GRADIENT";
		return NULL;
	} else if (scale > 0) {
		opt->done = false;
		opt->errmsg = NULL;

		vector_assign_copy(&opt->search_dir, grad0);
		vector_scale(&opt->search_dir, -1.0 / scale);
		double g0 = vector_dot(&opt->grad0, &opt->search_dir);
		assert(g0 < 0);
		return linesearch(opt, step0, f0, g0);
	} else {
		opt->done = true;
		opt->errmsg = NULL;
		return NULL;
	}
}

static bool converged(const struct bfgs *opt, double f,
		      const struct vector *grad)
{
	double abstol = opt->ctrl.abstol;
	double reltol = opt->ctrl.reltol;

	if (vector_norm(grad) < reltol * (vector_norm(&opt->x) + abstol)) {
		return true;
	}

	return false;
}

const struct vector *bfgs_advance(struct bfgs *opt, double f,
				  const struct vector *grad)
{
	assert(isfinite(f));
	assert(isfinite(vector_norm(grad)));

	double g = vector_dot(grad, &opt->search_dir);
	enum linesearch_task task = linesearch_advance(&opt->ls, f, g);
	bool wolfe = linesearch_sdec(&opt->ls) && linesearch_curv(&opt->ls);
	
	switch (task) {
	case LINESEARCH_CONV:
		break;

	case LINESEARCH_STEP:
		opt->ls_it++;

		if (opt->ls_it < opt->ctrl.ls_maxit) {
			vector_assign_copy(&opt->x, &opt->x0);
			vector_axpy(linesearch_step(&opt->ls), &opt->search_dir, &opt->x);
			return &opt->x;
		} else if (wolfe) {
			break;
		} else {
			opt->done = true;
			opt->errmsg = "linesearch failed to converge";
			return NULL;
		}
	default:
		if (wolfe) {
			break;
		} else {
			opt->done = true;
			opt->errmsg = linesearch_warnmsg(task);
			return NULL;
		}
	}

	update_hess(opt, grad);
	update_searchdir(opt, grad);

	// test for convergence
	opt->done = converged(opt, f, grad);

	vector_assign_copy(&opt->x0, &opt->x);
	vector_assign_copy(&opt->grad0, grad);
	double step0 = 1.0;
	double f0 = f;	
	double g0 = vector_dot(grad, &opt->search_dir);
	assert(g0 < 0);

	if (opt->done) {
		opt->errmsg = NULL;
		return NULL;
	} else {
		return linesearch(opt, step0, f0, g0);
	}
}

bool bfgs_converged(const struct bfgs *opt)
{
	return opt->done && !opt->errmsg;
}

const char *bfgs_error(const struct bfgs *opt)
{
	if (!opt->done) {
		return "BFGS FAILED TO CONVERGE";
	} else {
		return opt->errmsg;
	}
}
