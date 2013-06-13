#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "xalloc.h"

#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"

#include "linesearch.h"
#include "recv_fit.h"

/*
size_t constr_add_identify_recv_fit(struct constr *c, struct frame *f,
				    const struct messages *xmsgs,
				    const struct messages *ymsgs)
{
	struct recv_model model;
	struct recv_fit_eval eval;

	recv_model_init(&model, f, NULL);
	eval_init(&eval, &model, c);
	eval_set(&eval, NULL, xmsgs, ymsgs, f, &model, 2);

	size_t n = recv_coefs_dim(f);
	double *imatp = xcalloc(n * (n + 1) / 2, sizeof(double));

	recv_loglik_axpy_imat(1.0, &eval.loglik, imatp);

	size_t nadd = constr_add_identify(c, imatp, RECV_LOGLIK_UPLO);

	free(imatp);
	eval_deinit(&eval);
	recv_model_deinit(&model);

	return nadd;
}
*/
