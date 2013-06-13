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

	fit->imat = xmalloc(dim * (dim + 1) / 2 * sizeof(double));
	fit->uplo = RECV_LOGLIK_UPLO;

	memset(fit->imat, 0, dim * (dim + 1) / 2 * sizeof(double));
	recv_loglik_add_all(&fit->loglik, fit->msgs, fit->nmsg);
	recv_loglik_axpy_last_imat(1.0, &fit->loglik, fit->imat);

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
	recv_loglik_deinit(&fit->loglik);
	recv_model_deinit(&fit->model);
}
