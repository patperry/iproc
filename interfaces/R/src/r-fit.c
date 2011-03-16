
#include <stdio.h>
#include <R_ext/Rdynload.h>
#include "fit.h"
#include "r-messages.h"
#include "r-model.h"
#include "r-fit.h"


static R_CallMethodDef callMethods[] = {
    { "Riproc_fit", (DL_FUNC) &Riproc_fit, 3 },
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
            SEXP Rpenalty)
{
    iproc_model *model0 = Riproc_to_model(Rmodel0);
    iproc_messages *messages = Riproc_to_messages(Rmessages);
    double penalty = REAL(Rpenalty)[0];

    if(iproc_messages_max_from(messages) >= iproc_model_nsender(model0)) {
        error("message from id outside sender range");
    } else if (iproc_messages_max_to(messages) >= iproc_model_nreceiver(model0)) {
        error("message to id outside receiver range");
    } else if (!(penalty >= 0.0) && isfinite(penalty)
               && GET_LENGTH(Rpenalty) == 1) {
        error("penalty must be a positive finite scalar");
    }
    
    iproc_fit *fit = iproc_fit_new(model0, messages, penalty);
    int it = 0;
    
    do {
        it++;
        iproc_fit_step(fit);
        printf(".");
        fflush(stdout);
    } while (!iproc_fit_converged(fit, IPROC_FIT_ABSTOL, IPROC_FIT_RELTOL));
    
    
    SEXP Rmodel;
    PROTECT(Rmodel = Riproc_from_model(fit->model));
    iproc_fit_free(fit);

    UNPROTECT(1);
    return Rmodel;
}
