#ifndef _RIPROC_VARS_H
#define _RIPROC_VARS_H

#include <Rinternals.h>
#include "vars.h"

/* Call once to initialize library */
void         Riproc_vars_init ();

/* External functions */
SEXP         Riproc_vars_new       (SEXP Rsenders,
                                    SEXP Rreceivers);
SEXP         Riproc_vars_dim       (SEXP Rvars);
SEXP         Riproc_vars_nsender   (SEXP Rvars);
SEXP         Riproc_vars_nreceiver (SEXP Rvars);
SEXP         Riproc_vars_senders   (SEXP Rvars);
SEXP         Riproc_vars_receivers (SEXP Rvars);

/* Internal use only */
iproc_vars * Riproc_to_vars        (SEXP Rvars);
SEXP         Riproc_from_vars      (iproc_vars *vars);

#endif /* _RIPROC_VARS_H */
