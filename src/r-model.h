#ifndef _RIPROC_MODEL_H
#define _RIPROC_MODEL_H

#include <Rinternals.h>
#include "model.h"

/* Call once to initialize library */
void          Riproc_model_init ();

/* External functions */
SEXP          Riproc_model_new       (SEXP Rvars,
                                      SEXP Rcoefs,
                                      SEXP Rhas_loops);
SEXP          Riproc_model_vars      (SEXP Rmodel);
SEXP          Riproc_model_coefs     (SEXP Rmodel);
SEXP          Riproc_model_has_loops (SEXP Rmodel);

SEXP          Riproc_model_dim       (SEXP Rmodel);
SEXP          Riproc_model_nsender   (SEXP Rmodel);
SEXP          Riproc_model_nreceiver (SEXP Rmodel);

/* Internal use only */
iproc_model * Riproc_to_model        (SEXP Rmodel);
SEXP          Riproc_from_model      (iproc_model *model);

#endif /* _RIPROC_MODEL_H */
