#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "xalloc.h"

#include "recv_fit.h"


void recv_fit_init(struct recv_fit *fit, struct design *r, struct design2 *d,
		   int exclude_loops, const struct message *msgs, size_t nmsg,
		   const struct constr *c, const struct newton_ctrl *ctrl)
{
	struct newton_ctrl myctrl = NEWTON_CTRL0;
	size_t dim = design_dim(r) + design2_dim(d);

	assert(!c || constr_dim(c) == dim);

	fit->recv = r;
	fit->dyad = d;
	fit->msgs = msgs;
	fit->nmsg = nmsg;

	recv_model_init(&fit->model, NULL, fit->recv, fit->dyad);
	assert(dim == recv_model_dim(&fit->model));
	recv_model_set_exclude_loops(&fit->model, exclude_loops);

	recv_loglik_init(&fit->loglik, &fit->model);

	fit->nscore = xmalloc(dim * sizeof(double));
	fit->imat = xmalloc(dim * (dim + 1) / 2 * sizeof(double));
	fit->uplo = RECV_LOGLIK_UPLO;

	recv_loglik_add_all(&fit->loglik, fit->msgs, fit->nmsg);
	fit->df = recv_loglik_count(&fit->loglik);
	fit->dev0 = recv_loglik_dev(&fit->loglik);

	memset(fit->imat, 0, dim * (dim + 1) / 2 * sizeof(double));
	recv_loglik_axpy_imat(1.0, &fit->loglik, fit->imat);

	if (c) {
		constr_init_copy(&fit->constr, c);
	} else {
		constr_init(&fit->constr, dim);
	}
	fit->ncextra = constr_add_identify(&fit->constr, fit->imat, fit->uplo);

	myctrl.gtol *= MAX(1, fit->df * sqrt(dim));

	newton_init(&fit->opt, dim, &fit->constr, ctrl ? ctrl : &myctrl);
}


void recv_fit_deinit(struct recv_fit *fit)
{
	newton_deinit(&fit->opt);
	constr_deinit(&fit->constr);
	free(fit->imat);
	free(fit->nscore);
	recv_loglik_deinit(&fit->loglik);
	recv_model_deinit(&fit->model);
}


enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct recv_params *p,
				  const double *duals)
{
	const struct design *r = fit->recv;
	const struct design2 *d = fit->dyad;
	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_tvar_dim(r);
	size_t dimr = dimr0 + dimr1;
	size_t dimd0 = design2_trait_dim(d);
	size_t dimd1 = design2_tvar_dim(d);
	size_t dim = dimr + dimd0 + dimd1;

	assert(dim == recv_model_dim(&fit->model));

	double *x = xcalloc(dim, sizeof(double));
	if (p) {
		if (p->recv.traits)
			memcpy(x, p->recv.traits, dimr0 * sizeof(double));
		if (p->recv.tvars)
			memcpy(x + dimr0, p->recv.tvars, dimr1 * sizeof(double));
		if (p->dyad.traits)
			memcpy(x + dimr, p->dyad.traits, dimd0 * sizeof(double));
		if (p->dyad.tvars)
			memcpy(x + dimr + dimd0, p->dyad.tvars, dimd1 * sizeof(double));
	}

	recv_loglik_clear(&fit->loglik);
	recv_model_set_params(&fit->model, p);
	recv_model_set_moments(&fit->model, 2);
	recv_loglik_add_all(&fit->loglik, fit->msgs, fit->nmsg);
	fit->dev = recv_loglik_dev(&fit->loglik);

	struct recv_params nscore;
	nscore.recv.traits = fit->nscore;
	nscore.recv.tvars = fit->nscore + dimr0;
	nscore.dyad.traits = fit->nscore + dimr;
	nscore.dyad.tvars = fit->nscore + dimr + dimd0;

	memset(fit->nscore, 0, dim * sizeof(double));
	recv_loglik_axpy_score(-1.0, &fit->loglik, &nscore);

	memset(fit->imat, 0, dim * (dim + 1) / 2 * sizeof(double));
	recv_loglik_axpy_imat(1.0, &fit->loglik, fit->imat);

	fit->task = newton_start(&fit->opt, x, fit->dev, fit->nscore, duals);
	free(x);

	if (fit->task == NEWTON_ERR_DOM) {
		return RECV_FIT_ERR_DOM;
	} else if (fit->task == NEWTON_CONV) {
		return RECV_FIT_CONV;
	}
	assert(fit->task == NEWTON_HESS);

	fit->task = newton_set_hess(&fit->opt, fit->imat, fit->uplo);
	if (fit->task == NEWTON_ERR_HESS) {
		return RECV_FIT_ERR_IMAT;
	}
	assert(fit->task == NEWTON_STEP);
	while (fit->task == NEWTON_STEP) {
		const double *s = newton_search(&fit->opt);
		size_t ismax = blas_idamax(dim, s, 1) - 1;
		double smax = fabs(s[ismax]);
		double stpmax = 4.0 / MAX(1.0, smax);
		double stp = newton_step_size(&fit->opt);

		if (stp <= stpmax) {
			break;
		}

		fit->task = newton_step(&fit->opt, NAN, NULL);
	}

	return RECV_FIT_STEP;
}


enum recv_fit_task recv_fit_advance(struct recv_fit *fit)
{
	assert(fit->task == NEWTON_STEP || fit->task == NEWTON_HESS);

	const struct design *r = fit->recv;
	const struct design2 *d = fit->dyad;
	size_t dimr0 = design_trait_dim(r);
	size_t dimr1 = design_tvar_dim(r);
	size_t dimr = dimr0 + dimr1;
	size_t dimd0 = design2_trait_dim(d);
	size_t dimd1 = design2_tvar_dim(d);
	size_t dim = dimr + dimd0 + dimd1;
	enum recv_fit_task task;

	if (fit->task == NEWTON_HESS) {
		recv_model_set_moments(&fit->model, 2);
	} else {
		double *x = (double *)newton_next(&fit->opt);
		struct recv_params p;

		p.recv.traits = x;
		p.recv.tvars = x + dimr0;
		p.dyad.traits = x + dimr;
		p.dyad.tvars = x + dimr + dimd0;

		recv_model_set_moments(&fit->model, 1);
		recv_model_set_params(&fit->model, &p);
	}

	recv_loglik_clear(&fit->loglik);
	recv_loglik_add_all(&fit->loglik, fit->msgs, fit->nmsg);

	if (fit->task == NEWTON_STEP) {
		fit->dev = recv_loglik_dev(&fit->loglik);

		struct recv_params nscore;
		nscore.recv.traits = fit->nscore;
		nscore.recv.tvars = fit->nscore + dimr0;
		nscore.dyad.traits = fit->nscore + dimr;
		nscore.dyad.tvars = fit->nscore + dimr + dimd0;

		memset(fit->nscore, 0, dim * sizeof(double));
		recv_loglik_axpy_score(-1.0, &fit->loglik, &nscore);

		fit->task = newton_step(&fit->opt, fit->dev, fit->nscore);
	} else {
		memset(fit->imat, 0, dim * (dim + 1) / 2 * sizeof(double));
		recv_loglik_axpy_imat(1.0, &fit->loglik, fit->imat);

		fit->task = newton_set_hess(&fit->opt, fit->imat, fit->uplo);

		while (fit->task == NEWTON_STEP) {
			const double *s = newton_search(&fit->opt);
			size_t ismax = blas_idamax(dim, s, 1) - 1;
			double smax = fabs(s[ismax]);
			double stpmax = 4.0 / MAX(1.0, smax);
			double stp = newton_step_size(&fit->opt);

			if (stp <= stpmax) {
				break;
			}

			fit->task = newton_step(&fit->opt, NAN, NULL);
		}
	}

	switch (fit->task) {
	case NEWTON_CONV:
		task = RECV_FIT_CONV;
		break;
	case NEWTON_STEP:
	case NEWTON_HESS:
		task = RECV_FIT_STEP;
		break;
	case NEWTON_ERR_LNSRCH:
		task = RECV_FIT_ERR_LNSRCH;
		break;
	case NEWTON_ERR_XTOL:
		task = RECV_FIT_ERR_XTOL;
		break;
	case NEWTON_ERR_HESS:
		task = RECV_FIT_ERR_IMAT;
		break;
	case NEWTON_ERR_DOM:
		assert(0);
		task = RECV_FIT_ERR_DOM;
		break;
	}

	return task;
}

const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit)
{
	return &fit->loglik;
}

const struct constr *recv_fit_constr(const struct recv_fit *fit)
{
	return &fit->constr;
}

const double *recv_fit_duals(const struct recv_fit *fit)
{
	return newton_duals(&fit->opt);
}

double recv_fit_dev(const struct recv_fit *fit)
{
	return newton_val(&fit->opt);
}

void recv_fit_get_imat(const struct recv_fit *fit, const double **imat, enum blas_uplo *uplo)
{
	*imat = fit->imat;
	*uplo = fit->uplo;
}

double recv_fit_dev0(const struct recv_fit *fit)
{
	return fit->dev0;
}

double recv_fit_score_norm(const struct recv_fit *fit)
{
	return newton_grad_norm(&fit->opt);
}

double recv_fit_step_size(const struct recv_fit *fit)
{
	return newton_step_size(&fit->opt);
}

const char *recv_fit_errmsg(enum recv_fit_task task)
{
	switch (task) {
		case RECV_FIT_STEP:
			return "optimization in progress";
		case RECV_FIT_ERR_LNSRCH:
			return "linesearch failed";
		case RECV_FIT_ERR_XTOL:
			return "step size is less than tolerance";
		case RECV_FIT_ERR_IMAT:
			return "observed information matrix is singular";
		case RECV_FIT_ERR_DOM:
			return "initial point is not in domain";
		case RECV_FIT_CONV:
			return NULL;
	}

	assert(0);
	return NULL;
}

