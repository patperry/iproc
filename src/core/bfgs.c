#include "port.h"
#include <assert.h>
#include <math.h>
#include "bfgs.h"

static const struct vector *linesearch(struct bfgs *opt, double step0)
{
	opt->ls_it = 0;
	opt->step = step0;
	vector_assign_copy(&opt->x, &opt->x0);
	vector_axpy(opt->step, &opt->search_dir, &opt->x);

	linesearch_start(&opt->ls, opt->f0, opt->g0, &opt->ctrl.ls);
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
	vector_scale(s, opt->step);

	/* compute y */
	vector_assign_copy(y, grad);
	vector_axpy(-1.0, &opt->grad0, y);

	double s_y = vector_dot(s, y);

	/* update H */
	if (opt->first_step) {
		double scale = 1.0;

		if (s_y > 0) {
			double s_s = vector_dot(s, s);
			assert(s_s > 0);
			scale = s_y / s_s;
		}

		ssize_t i, n = vector_dim(y);

		matrix_fill(H, 0.0);
		for (i = 0; i < n; i++) {
			matrix_set_item(H, i, i, scale);
		}
		opt->first_step = false;
	} else {
		assert(s_y > 0);
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
	opt->f0 = f0;
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
		opt->g0 = vector_dot(&opt->grad0, &opt->search_dir);
		assert(opt->g0 < 0);
		return linesearch(opt, step0);
	} else {
		opt->done = true;
		opt->errmsg = NULL;

		opt->g0 = 0;
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

	/* linesearch in progress */
	if ((opt->step = linesearch_advance(&opt->ls, opt->step, f, g))) {
		opt->ls_it++;

		if (opt->ls_it == opt->ctrl.ls_maxit) {
			opt->done = true;
			opt->errmsg = "LINESEARCH FAILED TO CONVERGE";
			return NULL;
		}
		vector_assign_copy(&opt->x, &opt->x0);
		vector_axpy(opt->step, &opt->search_dir, &opt->x);
		return &opt->x;
	} else if (!linesearch_converged(&opt->ls)) {
		opt->done = true;
		opt->errmsg = linesearch_error(&opt->ls);
		return NULL;
	}

	opt->step = linesearch_step(&opt->ls);
	assert(opt->step > 0);
	update_hess(opt, grad);
	update_searchdir(opt, grad);

	// test for convergence
	opt->done = converged(opt, f, grad);

	vector_assign_copy(&opt->x0, &opt->x);
	opt->f0 = f;
	vector_assign_copy(&opt->grad0, grad);
	opt->g0 = vector_dot(grad, &opt->search_dir);
	assert(opt->g0 < 0);

	if (opt->done) {
		opt->errmsg = NULL;
		return NULL;
	} else {
		return linesearch(opt, 1.0);
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
