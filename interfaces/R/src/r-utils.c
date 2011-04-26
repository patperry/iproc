#include "port.h"
#include <assert.h>
#include <limits.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "r-utils.h"

static R_CallMethodDef callMethods[] = {
	{"Riproc_hash_numeric", (DL_FUNC) & Riproc_hash_numeric, 1},
	{NULL, NULL, 0}
};

void Riproc_utils_init(DllInfo * info)
{
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

SEXP Riproc_hash_numeric(SEXP x)
{
	struct vector view = Riproc_vector_view_sexp(x);
	size_t hash = vector_hash(&view);
	int int_hash = (int)(hash % INT_MAX);
	return ScalarInteger(int_hash);
}

/* Adapted from Luke Tierney's Regular expressions package
 *   http://www.stat.uiowa.edu/~luke/R/regexp.html
 */
Rboolean Riproc_sexpisptr(SEXP s, SEXP tag)
{
	if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != tag)
		return FALSE;

	return TRUE;
}

void *Riproc_sexp2ptr(SEXP s, Rboolean null_ok, SEXP tag, char *type)
{
	void *p;
	if (!Riproc_sexpisptr(s, tag))
		error("bad %s pointer", type);
	p = R_ExternalPtrAddr(s);
	if (!null_ok && p == NULL)
		error("null %s pointer", type);
	return p;
}

SEXP Riproc_vector_new_copy(struct vector * vector)
{
	int n = (int)vector_dim(vector);
	SEXP Rvector;

	PROTECT(Rvector = NEW_NUMERIC(n));
	struct vector view = Riproc_vector_view_sexp(Rvector);
	vector_assign_copy(&view, vector);
	UNPROTECT(1);
	return Rvector;
}

SEXP Riproc_matrix_new_copy(struct matrix * matrix)
{
	int nrow = (int)matrix_nrow(matrix);
	int ncol = (int)matrix_ncol(matrix);
	SEXP Rmatrix;

	PROTECT(Rmatrix = allocMatrix(REALSXP, nrow, ncol));
	struct matrix view = Riproc_matrix_view_sexp(Rmatrix);
	matrix_assign_copy(&view, matrix);
	UNPROTECT(1);
	return Rmatrix;
}

struct vector Riproc_vector_view_sexp(SEXP Rvector)
{
	int n = GET_LENGTH(Rvector);
	double *data = NUMERIC_POINTER(Rvector);
	struct vector view;
	vector_init_view(&view, data, n);
	return view;
}

struct matrix Riproc_matrix_view_sexp(SEXP Rmatrix)
{
	struct matrix view;

	if (isMatrix(Rmatrix)) {
		int nrow = INTEGER(GET_DIM(Rmatrix))[0];
		int ncol = INTEGER(GET_DIM(Rmatrix))[1];
		double *data = NUMERIC_POINTER(Rmatrix);
		matrix_init_view(&view, data, nrow, ncol);
	} else {
		int nrow = GET_LENGTH(Rmatrix);
		int ncol = 1;
		double *data = NUMERIC_POINTER(Rmatrix);
		matrix_init_view(&view, data, nrow, ncol);
	}

	return view;
}
