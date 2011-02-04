#ifndef _RIPROC_FRAME_H
#define _RIPROC_FRAME_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "frame.h"
#include "vrecip.h"

/* Call once to initialize library */
void         Riproc_frame_init      (DllInfo *info);

/* External functions */
SEXP         Riproc_frame_new       (SEXP Rsenders,
                                     SEXP Rreceivers,
                                     SEXP Rreceive_intervals);
SEXP         Riproc_frame_dim       (SEXP Rframe);
SEXP         Riproc_frame_nsender   (SEXP Rframe);
SEXP         Riproc_frame_nreceiver (SEXP Rframe);
SEXP         Riproc_frame_senders   (SEXP Rframe);
SEXP         Riproc_frame_receivers (SEXP Rframe);

SEXP         Riproc_frame_mul       (SEXP Rframe,
                                     SEXP Rmatrix,
                                     SEXP Rsender,
                                     SEXP Rcursor);
SEXP         Riproc_frame_tmul      (SEXP Rframe,
                                     SEXP Rmatrix,
                                     SEXP Rsender,
                                     SEXP Rcursor);


/* Internal use only */
iproc_frame * Riproc_to_frame       (SEXP Rframe);
SEXP          Riproc_from_frame     (iproc_frame *frame);

typedef struct _Riproc_frame_udata Riproc_frame_udata;

struct _Riproc_frame_udata {
    iproc_vrecip *recip;
};

Riproc_frame_udata * Riproc_frame_udata_new             (SEXP               Rreceive_intervals);
void                 Riproc_frame_udata_free             (void              *frame_udata);
int64_t              Riproc_frame_udata_dim              (Riproc_frame_udata *udata);
void                 Riproc_frame_udata_get_sender_frame (iproc_frame_ctx    *ctx);

#endif /* _RIPROC_FRAME_H */
