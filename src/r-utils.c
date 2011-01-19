
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "r-utils.h"

/* Adapted from Luke Tierney's Regular expressions package
 *   http://www.stat.uiowa.edu/~luke/R/regexp.html
 */
Rboolean
Riproc_sexpisptr (SEXP s,
                  SEXP tag)
{
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != tag)
        return FALSE;

    return TRUE;
}


void *
Riproc_sexp2ptr (SEXP     s,
                 Rboolean null_ok,
                 SEXP     tag,
                 char    *type)
{
    void *p;
    if (!Riproc_sexpisptr(s, tag))
        error("bad %s pointer", type);
    p = R_ExternalPtrAddr(s);
    if (!null_ok && p == NULL)
        error("null %s pointer", type);
    return p;
}

SEXP
Riproc_vector_new_copy (iproc_vector *vector)
{
    int64_t n = iproc_vector_dim(vector);
    SEXP Rvector;

    PROTECT(Rvector = NEW_NUMERIC(n));
    iproc_vector_view view = Riproc_vector_view_sexp(Rvector);
    iproc_vector_copy(&view.vector, vector);
    UNPROTECT(1);
    return Rvector;
}

SEXP
Riproc_matrix_new_copy (iproc_matrix *matrix)
{
    int64_t nrow = iproc_matrix_nrow(matrix);
    int64_t ncol = iproc_matrix_ncol(matrix);
    SEXP Rmatrix;

    PROTECT(Rmatrix = allocMatrix(REALSXP, nrow, ncol));
    iproc_matrix_view view = Riproc_matrix_view_sexp(Rmatrix);
    iproc_matrix_copy(&view.matrix, matrix);
    UNPROTECT(1);
    return Rmatrix;
}

iproc_vector_view
Riproc_vector_view_sexp (SEXP Rvector)
{
    int64_t n = GET_LENGTH(Rvector);
    double *data = NUMERIC_POINTER(Rvector);
    return iproc_vector_view_array(data, n);
}

iproc_matrix_view
Riproc_matrix_view_sexp (SEXP Rmatrix)
{
    if (isMatrix(Rmatrix)) {
        int64_t nrow = INTEGER(GET_DIM(Rmatrix))[0];
        int64_t ncol = INTEGER(GET_DIM(Rmatrix))[1];
        double *data = NUMERIC_POINTER(Rmatrix);
        return iproc_matrix_view_array(data, nrow, ncol);
    } else {
        int64_t nrow = GET_LENGTH(Rmatrix);
        int64_t ncol = 1;
        double *data = NUMERIC_POINTER(Rmatrix);
        return iproc_matrix_view_array(data, nrow, ncol);
    }
}

