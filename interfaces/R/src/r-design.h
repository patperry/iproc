#ifndef _RIPROC_DESIGN_H
#define _RIPROC_DESIGN_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "design.h"

/* Call once to initialize library */
void Riproc_design_init(DllInfo * info);

/* External functions */
SEXP Riproc_design_new(SEXP Rsenders,
		       SEXP Rreceivers,
		       SEXP Rreceiver_effects, SEXP Rrecip_intervals);
SEXP Riproc_design_dim(SEXP Rdesign);
SEXP Riproc_design_nsender(SEXP Rdesign);
SEXP Riproc_design_nreceiver(SEXP Rdesign);
SEXP Riproc_design_senders(SEXP Rdesign);
SEXP Riproc_design_receivers(SEXP Rdesign);

/* Internal use only */
const struct design *Riproc_to_design(SEXP Rdesign);
SEXP Riproc_from_design(const struct design *design);

#endif /* _RIPROC_DESIGN_H */
