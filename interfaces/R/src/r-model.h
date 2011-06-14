#ifndef _RIPROC_MODEL_H
#define _RIPROC_MODEL_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "model.h"

/* Call once to initialize library */
void Riproc_model_init(DllInfo * info);

/* External functions */
SEXP Riproc_model_new(SEXP Rdesign, SEXP Rcoefs);
SEXP Riproc_model_design(SEXP Rmodel);
SEXP Riproc_model_coefs(SEXP Rmodel);

SEXP Riproc_model_dim(SEXP Rmodel);
SEXP Riproc_model_nsender(SEXP Rmodel);
SEXP Riproc_model_nreceiver(SEXP Rmodel);

// SEXP Riproc_model_log_probs(SEXP Rmodel, SEXP Risend, SEXP Rcursor);

/* Internal use only */
iproc_model *Riproc_to_model(SEXP Rmodel);
SEXP Riproc_from_model(iproc_model * model);

#endif /* _RIPROC_MODEL_H */
