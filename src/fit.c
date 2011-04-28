#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fit.h"

int dcsrch_(double *stp, double *f, double *g,
	    double *ftol, double *gtol, double *xtol, char *task,
	    double *stpmin, double *stpmax, f77int * isave, double *dsave,
	    f77int task_len);

static void
eval_objective(iproc_loglik * loglik,
	       double penalty, double *valuep, struct vector *grad)
{
	iproc_model *model = loglik->model;
	struct vector *coefs = model->coefs;

	double n = loglik->nrecv;
	double ll_value = iproc_loglik_value(loglik);
	struct vector *ll_grad = iproc_loglik_grad(loglik);

	double norm = vector_norm(coefs);
	double norm2 = norm * norm;
	double value = -ll_value / n + 0.5 * (penalty / n) * norm2;
	*valuep = value;

	vector_assign_copy(grad, coefs);
	vector_scale(grad, penalty / n);
	vector_axpy(-1.0 / n, ll_grad, grad);
}

static bool iproc_fit_init(iproc_fit * fit)
{
	ssize_t dim = iproc_model_dim(fit->model);

	fit->inv_hess = NULL;
	fit->x0 = vector_alloc(dim);
	fit->x = vector_alloc(dim);
	fit->grad0 = vector_alloc(dim);
	fit->grad = vector_alloc(dim);
	fit->search_dir = vector_alloc(dim);
	fit->loglik = iproc_loglik_new(fit->model, fit->messages);

	if (!
	    (fit->x0 && fit->x && fit->grad0 && fit->grad && fit->search_dir
	     && fit->loglik))
		return false;

	eval_objective(fit->loglik, fit->penalty, &fit->value, fit->grad);
	vector_assign_copy(fit->x, fit->model->coefs);
	fit->value0 = NAN;
	fit->step = 1.0;

	return true;
}

iproc_fit *iproc_fit_new(iproc_model * model0,
			 struct messages * messages, double penalty)
{
	assert(model0);
	assert(messages);
	assert(penalty >= 0.0);
	assert(iproc_messages_max_from(messages) < iproc_model_nsender(model0));
	assert(iproc_messages_max_to(messages) < iproc_model_nreceiver(model0));

	iproc_fit *fit = malloc(sizeof(*fit));
	if (!fit)
		return NULL;

	fit->model = iproc_model_ref(model0);
	fit->messages = messages_ref(messages);
	fit->penalty = penalty;
	fit->loglik = NULL;

	if (!iproc_fit_init(fit)) {
		iproc_fit_free(fit);
		fit = NULL;
	}

	return fit;
}

void iproc_fit_free(iproc_fit * fit)
{
	if (fit) {
		iproc_loglik_unref(fit->loglik);
		vector_free(fit->search_dir);
		vector_free(fit->grad);
		vector_free(fit->grad0);
		vector_free(fit->x);
		vector_free(fit->x0);
		matrix_free(fit->inv_hess);
		messages_free(fit->messages);
		iproc_model_unref(fit->model);
		free(fit);
	}
}

static void linesearch(iproc_fit * fit)
{
	assert(fit);

	iproc_model *model = fit->model;
	struct messages *messages = fit->messages;
	double penalty = fit->penalty;
	iproc_design *design = iproc_design_ref(model->design);
	int has_loops = model->has_loops;
	iproc_loglik *loglik = fit->loglik;

	struct vector *search_dir = fit->search_dir;
	struct vector *grad0 = fit->grad;	/* swap grad0 and grad */
	struct vector *grad = fit->grad0;
	struct vector *x0 = fit->x;	/* swap x0 and x */
	struct vector *x = fit->x0;
	double value0 = fit->value;
	double value;

	double f = value0;
	double g = vector_dot(grad0, search_dir);
	double ftol = 1e-4;
	double gtol = 0.9;
	double xtol = 0.1;
	f77int task_len = 60;
	char task[task_len + 1];
	double stpmin = 0;
	double stpmax = -f / (gtol * g);
	double stp = g == 0 ? 1.0 : MIN(1.0, stpmax);
	f77int isave[2];
	double dsave[13];

	if (f == 0 || g == 0)
		goto cleanup;

	strcpy(task, "START");
	dcsrch_(&stp, &f, &g, &ftol, &gtol, &xtol, task, &stpmin, &stpmax,
		isave, dsave, task_len);

	while (strncmp(task, "FG", strlen("FG")) == 0) {
		/* Take a step and create a new model */
		vector_assign_copy(x, x0);
		vector_axpy(stp, search_dir, x);

		iproc_model_unref(model);
		model = iproc_model_new(design, x, has_loops);

		/* Update the loglik, value, and gradient */
		iproc_loglik_unref(loglik);
		loglik = iproc_loglik_new(model, messages);
		eval_objective(loglik, penalty, &value, grad);

		f = value;
		g = vector_dot(grad, search_dir);

		dcsrch_(&stp, &f, &g, &ftol, &gtol, &xtol, task, &stpmin,
			&stpmax, isave, dsave, task_len);
	}

	const char *xtol_msg = "WARNING: XTOL TEST SATISFIED";
	if (!(strncmp(task, "CONV", strlen("CONV")) == 0
	      || strncmp(task, xtol_msg, strlen(xtol_msg)) == 0)) {
		printf("%s\n", task);
	}

cleanup:
	iproc_design_unref(design);

	fit->step = stp;
	fit->x = x;
	fit->x0 = x0;
	fit->value = value;
	fit->value0 = value0;
	fit->grad = grad;
	fit->grad0 = grad0;
	fit->model = model;
	fit->loglik = loglik;
}

static void update_hess(iproc_fit * fit)
{
	struct vector *s = vector_alloc_copy(fit->search_dir);
	vector_scale(s, fit->step);

	struct vector *y = vector_alloc_copy(fit->grad);
	vector_axpy(-1.0, fit->grad0, y);

	double s_y = vector_dot(s, y);

	struct matrix *H = fit->inv_hess;

	if (H == NULL) {
		double scale = 1.0;

		if (s_y > 0) {
			double s_s = vector_dot(s, s);
			scale = s_y / s_s;
		}

		ssize_t i, n = vector_dim(y);
		H = matrix_alloc(n, n);
		matrix_fill(H, 0.0);

		for (i = 0; i < n; i++) {
			matrix_set(H, i, i, scale);
		}

		fit->inv_hess = H;
	} else {
		struct vector *H_y = vector_alloc(vector_dim(y));
		matrix_mul(1.0, TRANS_NOTRANS, H, y, 0.0, H_y);

		double y_H_y = vector_dot(H_y, y);
		double scale1 = (1.0 + (y_H_y / s_y)) / s_y;
		double rho = 1.0 / s_y;

		matrix_update1(H, scale1, s, s);
		matrix_update1(H, -rho, H_y, s);
		matrix_update1(H, -rho, s, H_y);

		vector_free(H_y);
	}

	vector_free(y);
	vector_free(s);
}

static void update_searchdir(iproc_fit * fit)
{
	assert(fit);

	struct vector *s = fit->search_dir;

	if (fit->inv_hess) {
		matrix_mul(-1.0, TRANS_NOTRANS, fit->inv_hess,
			   fit->grad, +0.0, s);
	} else {
		vector_assign_copy(s, fit->grad);
		double scale = vector_norm(s);
		ssize_t i, n = vector_dim(s);

		if (scale != 0) {
			for (i = 0; i < n; i++) {
				double s_i = *vector_at(s, i);
				*vector_at(s, i) = -s_i / scale;
			}
		}
	}
}

void iproc_fit_step(iproc_fit * fit)
{
	assert(fit);
	update_searchdir(fit);
	linesearch(fit);
	update_hess(fit);
}

bool iproc_fit_converged(iproc_fit * fit, double abs_tol, double rel_tol)
{
	double f = fit->value;
	double f0 = fit->value0;
	double step = fit->step;
	double g0 = vector_dot(fit->search_dir, fit->grad);
	bool converged = false;

	if (fabs(f - f0) <= abs_tol && step * fabs(g0) <= abs_tol) {
		/* abs_tol test satisfied */
		converged = true;
	} else if (fabs(f - f0) <= rel_tol * fabs(f0)
		   && step * fabs(g0) <= rel_tol * fabs(f0)) {
		/* rel_tol test satisfied */
		converged = true;
	}

	return converged;
}
