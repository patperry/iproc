#include "port.h"
#include <assert.h>

#include <stdio.h>
#include <R_ext/Print.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>

#include "fit.h"
#include "bfgs.h"
#include "r-messages.h"
#include "r-design.h"
#include "r-fit.h"

static R_CallMethodDef callMethods[] = {
	{"Riproc_fit", (DL_FUNC) & Riproc_fit, 8},
	{NULL, NULL, 0}
};

void Riproc_fit_init(DllInfo * info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

SEXP
Riproc_fit(SEXP Rdesign,
	   SEXP Rmessages,
	   SEXP Rpenalty,
	   SEXP Rreltol, SEXP Rabstol, SEXP Rmaxit, SEXP Rtrace, SEXP Rreport)
{
	struct design *design = Riproc_to_design(Rdesign);
	struct messages *messages = Riproc_to_messages(Rmessages);
	double penalty = REAL(Rpenalty)[0];
	double reltol = REAL(Rreltol)[0];
	double abstol = REAL(Rabstol)[0];
	int maxit = INTEGER_VALUE(Rmaxit);
	bool trace = LOGICAL_VALUE(Rtrace);
	int report = INTEGER_VALUE(Rreport);

	if (messages_max_from(messages) >= design_send_count(design)) {
		error("message 'from' id outside sender range");
	} else if (messages_max_to(messages) >= design_recv_count(design)) {
		error("message 'to' id outside receiver range");
	} else if (!(penalty >= 0.0 && isfinite(penalty))) {
		error("value of 'penalty' must be >= 0 and finite");
	} else if (!(reltol > 0)) {
		error("value of 'reltol' must be > 0");
	} else if (!(abstol > 0)) {
		error("value of 'abstol' must be > 0");
	} else if (maxit == NA_INTEGER || !(maxit > 0)) {
		error("value of 'maxit' must be > 0");
	} else if (trace == NA_LOGICAL) {
		error("'trace' value is NA");
	} else if (report == NA_INTEGER || !(report > 0)) {
		error("'report' value must be positive");
	}

	struct recv_fit fit;
	recv_fit_init(&fit, messages, design, NULL, penalty);
	int it = 0;

	do {
		R_CheckUserInterrupt();
		it++;
		recv_fit_step(&fit);

		if (trace && it % report == 0) {
			const char *msg = penalty == 0 ? "" : "(penalized) ";
			ssize_t n = recv_loglik_count(&fit.loglik);
			double dev = n * recv_loglik_avg_dev(&fit.loglik);
			double dec = -2 * bfgs_decrement(&fit.opt);
			Rprintf("iter %d deviance %s%.6f decrement %.6f\n", it,
				msg, dev, dec);
		}
	} while (it < maxit && !recv_fit_converged(&fit));

	if (it == maxit && !recv_fit_converged(&fit)) {
		warning("algorithm did not converge");
	}

	assert(0 && "return value not implemented");
	return NULL;
}
