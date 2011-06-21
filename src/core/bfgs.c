#include "port.h"
#include <assert.h>
#include <math.h>
#include "bfgs.h"


void bfgs_init(struct bfgs *opt, ssize_t n, const struct bfgs_ctrl *ctrl)
{
	assert(opt);
	assert(n >= 0);
	assert(ctrl);
	assert(bfgs_ctrl_valid(ctrl));
	
	opt->ctrl = *ctrl;
	matrix_init(&opt->inv_hess, n, n);
	vector_init(&opt->search_dir, n);
	vector_init(&opt->g0, n);
	vector_init(&opt->g1, n);
	vector_init(&opt->s, n);
	vector_init(&opt->y, n);
	vector_init(&opt->H_y, n);
	bfgs_clear(opt);
}

void bfgs_deinit(struct bfgs *opt)
{
	vector_deinit(&opt->H_y);	
	vector_deinit(&opt->y);
	vector_deinit(&opt->s);
	vector_deinit(&opt->g1);
	vector_deinit(&opt->g0);
	vector_deinit(&opt->search_dir);	
	matrix_deinit(&opt->inv_hess);
}

void bfgs_clear(struct bfgs *opt)
{
	opt->grad0 = &opt->g0;
	opt->grad = &opt->g1;
	opt->step = NAN;
	opt->started = false;	
	opt->swapped = false;
}

static void update_searchdir(struct bfgs *opt)
{
	assert(opt);
	
	struct vector *s = &opt->search_dir;
	
	if (opt->started) {
		matrix_mul(-1.0, TRANS_NOTRANS, &opt->inv_hess,
			   opt->grad, 0.0, s);
	} else {
		vector_assign_copy(s, opt->grad);
		double scale = vector_norm(s);
		
		if (scale != 0) {
			vector_scale(s, -1.0 / scale);
		}
	}
}


static void update_hess(struct bfgs *opt)
{
	struct matrix *H = &opt->inv_hess;
	struct vector *s = &opt->s;
	struct vector *y = &opt->y;

	/* compute s */
	vector_assign_copy(s, &opt->search_dir);
	vector_scale(s, opt->step);
	
	/* compute y */
	vector_assign_copy(y, opt->grad);
	vector_axpy(-1.0, opt->grad0, y);
	
	double s_y = vector_dot(s, y);
	
	/* update H */
	if (!opt->started) {
		double scale = 1.0;
		
		if (s_y > 0) {
			double s_s = vector_dot(s, s);
			scale = s_y / s_s;
		}
		
		ssize_t i, n = vector_dim(y);
		
		matrix_fill(H, 0.0);
		for (i = 0; i < n; i++) {
			matrix_set_item(H, i, i, scale);
		}
	} else {
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
