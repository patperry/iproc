#ifndef _RIPROC_UTILS_H
#define _RIPROC_UTILS_H

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#include "vector.h"
#include "matrix.h"


/* Call once to initialize library */
void              Riproc_utils_init       (DllInfo *info);

/* External functions */
SEXP              Riproc_hash_numeric     (SEXP     x);

/* Internal use only */
void *            Riproc_sexp2ptr         (SEXP     s,
                                           Rboolean null_ok,
                                           SEXP     tag,
                                           char    *type);
Rboolean          Riproc_sexpisptr        (SEXP     s,
                                           SEXP     tag);

SEXP              Riproc_vector_new_copy  (iproc_vector *vector);
SEXP              Riproc_matrix_new_copy  (iproc_matrix *matrix);

iproc_vector_view Riproc_vector_view_sexp (SEXP Rvector);
iproc_matrix_view Riproc_matrix_view_sexp (SEXP Rmatrix);

#endif /* _RIPROC_UTILS_H */
