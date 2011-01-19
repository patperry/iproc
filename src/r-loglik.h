#ifndef _RIPROC_LOGLIK_H
#define _RIPROC_LOGLIK_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "loglik.h"

/* Call once to initialize library */
void           Riproc_loglik_init  (DllInfo *info);

/* External functions */
SEXP           Riproc_loglik_new   (SEXP          Rmodel,
                                    SEXP          Rmessages);
SEXP           Riproc_loglik_value (SEXP          Rloglik);
SEXP           Riproc_loglik_grad  (SEXP          Rloglik);

/* Internal use only */
iproc_loglik * Riproc_to_loglik    (SEXP          Rloglik);
SEXP           Riproc_from_loglik  (iproc_loglik *loglik);

#endif /* _RIPROC_LOGLIK_H */
