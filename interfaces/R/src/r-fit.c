
#include <stdio.h>
#include <R_ext/Print.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Utils.h>

#include "fit.h"
#include "r-messages.h"
#include "r-model.h"
#include "r-fit.h"


static R_CallMethodDef callMethods[] = {
    { "Riproc_fit", (DL_FUNC) &Riproc_fit, 8 },
    { NULL,         NULL,                  0 }
};


void
Riproc_fit_init (DllInfo *info)
{
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

SEXP
Riproc_fit (SEXP Rmodel0,
            SEXP Rmessages,
            SEXP Rpenalty,
            SEXP Rreltol,
            SEXP Rabstol,
            SEXP Rmaxit,
            SEXP Rtrace,
            SEXP Rreport)
{
    iproc_model *model0 = Riproc_to_model(Rmodel0);
    iproc_messages *messages = Riproc_to_messages(Rmessages);
    double penalty = REAL(Rpenalty)[0];
    double reltol = REAL(Rreltol)[0];
    double abstol = REAL(Rabstol)[0];
    int maxit = INTEGER_VALUE(Rmaxit);
    bool trace = LOGICAL_VALUE(Rtrace);
    int report = INTEGER_VALUE(Rreport);
    

    if(iproc_messages_max_from(messages) >= iproc_model_nsender(model0)) {
        error("message 'from' id outside sender range");
    } else if (iproc_messages_max_to(messages) >= iproc_model_nreceiver(model0)) {
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
    
    iproc_fit *fit = iproc_fit_new(model0, messages, penalty);
    int it = 0;
    
    do {
        R_CheckUserInterrupt();
        it++;
        iproc_fit_step(fit);
        
        if (trace && it % report == 0) {
            const char * msg = penalty == 0 ? "" : "(penalized) ";
            int64_t n = fit->loglik->nrecv;
            double dev = 2 * fit->value * n;
            double dec = -1 * iproc_vector_dot(fit->search_dir, fit->grad) * n;
            Rprintf("iter %d deviance %s%.6f decrement %.6f\n",
                    it, msg, dev, dec);
        }
    } while (it < maxit
             && !iproc_fit_converged(fit, abstol, reltol));
    
    if (it == maxit && !iproc_fit_converged(fit, abstol, reltol)) {
        warning("algorithm did not converge");
    }
    
    
    SEXP Rmodel;
    PROTECT(Rmodel = Riproc_from_model(fit->model));
    iproc_fit_free(fit);

    UNPROTECT(1);
    return Rmodel;
}
