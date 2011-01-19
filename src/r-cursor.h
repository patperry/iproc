#ifndef _RIPROC_CURSOR_H
#define _RIPROC_CURSOR_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "cursor.h"


/* Call once to initialize library */
void           Riproc_cursor_init    (DllInfo      *info);

/* External functions */
SEXP           Riproc_cursor_new     (SEXP          Rmessages);
SEXP           Riproc_cursor_advance (SEXP          Rcursor);
SEXP           Riproc_cursor_reset   (SEXP          Rcursor);

SEXP           Riproc_cursor_time    (SEXP          Rcursor);
SEXP           Riproc_cursor_nties   (SEXP          Rcursor);
SEXP           Riproc_cursor_from    (SEXP          Rcursor);
SEXP           Riproc_cursor_to      (SEXP          Rcursor);

/* Internal use only */
iproc_cursor * Riproc_to_cursor      (SEXP          Rcursor);
SEXP           Riproc_from_cursor    (iproc_cursor *cursor);


#endif /* _RIPROC_CURSOR_H */
