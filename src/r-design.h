#ifndef _RIPROC_DESIGN_H
#define _RIPROC_DESIGN_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "design.h"
#include "vrecip.h"

/* Call once to initialize library */
void         Riproc_design_init      (DllInfo *info);

/* External functions */
SEXP         Riproc_design_new       (SEXP Rsenders,
                                     SEXP Rreceivers,
                                     SEXP Rrecip_intervals);
SEXP         Riproc_design_dim       (SEXP Rdesign);
SEXP         Riproc_design_nsender   (SEXP Rdesign);
SEXP         Riproc_design_nreceiver (SEXP Rdesign);
SEXP         Riproc_design_senders   (SEXP Rdesign);
SEXP         Riproc_design_receivers (SEXP Rdesign);

SEXP         Riproc_design_mul       (SEXP Rdesign,
                                     SEXP Rmatrix,
                                     SEXP Rsender,
                                     SEXP Rcursor);
SEXP         Riproc_design_tmul      (SEXP Rdesign,
                                     SEXP Rmatrix,
                                     SEXP Rsender,
                                     SEXP Rcursor);


/* Internal use only */
iproc_design * Riproc_to_design       (SEXP Rdesign);
SEXP          Riproc_from_design     (iproc_design *design);

typedef struct _Riproc_design_udata Riproc_design_udata;

struct _Riproc_design_udata {
    iproc_vrecip *recip;
};

Riproc_design_udata * Riproc_design_udata_new             (SEXP               Rrecip_intervals);
void                 Riproc_design_udata_free             (void              *design_udata);
int64_t              Riproc_design_udata_dim              (Riproc_design_udata *udata);
void                 Riproc_design_udata_get_sdesign_vars (iproc_design_ctx    *ctx);

#endif /* _RIPROC_DESIGN_H */
