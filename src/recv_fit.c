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
		   const struct message *msgs, size_t nmsg,
		   const struct constr *c, const struct newton_ctrl *ctrl)
{
	size_t dim = design_dim(r) + design2_dim(d);

	assert(!c || constr_dim(c) == dim);

	fit->recv = r;
	fit->dyad = d;
	fit->msgs = msgs;
	fit->nmsg = nmsg;

	recv_model_init(&fit->model, NULL, fit->recv, fit->dyad);
	assert(dim == recv_model_dim(&fit->model));

	recv_loglik_init(&fit->loglik, &fit->model);

	fit->score = xmalloc(dim * sizeof(double));
	fit->imat = xmalloc(dim * (dim + 1) / 2 * sizeof(double));
	fit->uplo = RECV_LOGLIK_UPLO;

	memset(fit->imat, 0, dim * (dim + 1) / 2 * sizeof(double));
	recv_loglik_add_all(&fit->loglik, fit->msgs, fit->nmsg);
	recv_loglik_axpy_imat(1.0, &fit->loglik, fit->imat);

	if (c) {
		constr_init_copy(&fit->constr, c);
	} else {
		constr_init(&fit->constr, dim);
	}
	fit->ncextra = constr_add_identify(&fit->constr, fit->imat, fit->uplo);

	newton_init(&fit->opt, dim, &fit->constr, ctrl);
}


void recv_fit_deinit(struct recv_fit *fit)
{
	newton_deinit(&fit->opt);
	constr_deinit(&fit->constr);
	free(fit->imat);
	free(fit->score);
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
	enum newton_task ntask;

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

	struct recv_params score;
	score.recv.traits = fit->score;
	score.recv.tvars = fit->score + dimr0;
	score.dyad.traits = fit->score + dimr;
	score.dyad.tvars = fit->score + dimr + dimd0;

	memset(fit->score, 0, dim * sizeof(double));
	recv_loglik_axpy_score(1.0, &fit->loglik, &score);

	memset(fit->imat, 0, dim * (dim + 1) / 2 * sizeof(double));
	recv_loglik_axpy_imat(1.0, &fit->loglik, fit->imat);

	ntask = newton_start(&fit->opt, x, fit->dev, fit->score, duals);
	if (ntask == NEWTON_ERR_DOM) {
		return RECV_FIT_ERR_DOM;
	} else if (ntask == NEWTON_CONV) {
		return RECV_FIT_CONV;
	}
	assert(ntask == NEWTON_HESS);

	ntask = newton_set_hess(&fit->opt, fit->imat, fit->uplo);
	if (ntask == NEWTON_ERR_HESS) {
		return RECV_FIT_ERR_IMAT;
	}
	assert(ntask == NEWTON_STEP);

	return RECV_FIT_STEP;
}



