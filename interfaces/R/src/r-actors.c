#include "port.h"
#include <assert.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "actors.h"
#include "vector.h"
#include "matrix.h"
#include "r-utils.h"
#include "r-actors.h"

static SEXP Riproc_actors_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_actors_new", (DL_FUNC) & Riproc_actors_new, 1},
	{"Riproc_actors_size", (DL_FUNC) & Riproc_actors_size, 1},
	{"Riproc_actors_dim", (DL_FUNC) & Riproc_actors_dim, 1},
	{"Riproc_actors_traits", (DL_FUNC) & Riproc_actors_traits, 2},
	{"Riproc_actors_mul", (DL_FUNC) & Riproc_actors_mul, 2},
	{"Riproc_actors_tmul", (DL_FUNC) & Riproc_actors_tmul, 2},
	{NULL, NULL, 0}
};

void Riproc_actors_init(DllInfo * info)
{
	Riproc_actors_type_tag = install("Riproc_actors_type_tag");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void Riproc_actors_free(SEXP Ractors)
{
	struct actors *actors = Riproc_to_actors(Ractors);
	actors_free(actors);
}

SEXP Riproc_from_actors(struct actors *actors)
{
	SEXP Ractors, class;

	actors_ref(actors);

	/* store the actors pointer in an external pointer */
	PROTECT(Ractors =
		R_MakeExternalPtr(actors, Riproc_actors_type_tag, R_NilValue));
	R_RegisterCFinalizer(Ractors, Riproc_actors_free);

	/* set the class of the result */
	PROTECT(class = allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("actors"));
	classgets(Ractors, class);

	UNPROTECT(2);
	return Ractors;
}

struct actors *Riproc_to_actors(SEXP Ractors)
{
	Rboolean null_ok = TRUE;
	SEXP tag = Riproc_actors_type_tag;
	char *type = "actors";
	struct actors *actors = Riproc_sexp2ptr(Ractors, null_ok, tag, type);
	return actors;
}

SEXP Riproc_actors_new(SEXP Rtraits_t)
{
	struct matrix traits_t = Riproc_matrix_view_sexp(Rtraits_t);
	int n = (int)matrix_ncol(&traits_t);

	if (!(n > 0))
		error("must specify at least one actor");

	struct vector traits0 = matrix_col(&traits_t, 0);
	ssize_t dim = vector_dim(&traits0);
	struct actors *actors = actors_alloc(dim);
	SEXP Ractors;

	if (!actors)
		error("could not allocate new actors object");

	actors_init_matrix(actors, &traits_t, TRANS_TRANS);

	PROTECT(Ractors = Riproc_from_actors(actors));
	actors_free(actors);
	UNPROTECT(1);
	return Ractors;
}

SEXP Riproc_actors_size(SEXP Ractors)
{
	struct actors *actors = Riproc_to_actors(Ractors);
	int size = (int)actors_size(actors);
	return ScalarInteger(size);
}

SEXP Riproc_actors_dim(SEXP Ractors)
{
	struct actors *actors = Riproc_to_actors(Ractors);
	int dim = (int)actors_dim(actors);
	return ScalarInteger(dim);
}

SEXP Riproc_actors_traits(SEXP Ractors, SEXP Ractor_ids)
{
	struct actors *actors = Riproc_to_actors(Ractors);
	int size = (int)actors_size(actors);
	int dim = (int)actors_dim(actors);
	int i, n = GET_LENGTH(Ractor_ids);
	SEXP Rxt;

	PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
	struct matrix xt = Riproc_matrix_view_sexp(Rxt);
	struct vector dst;
	const struct vector *src;

	for (i = 0; i < n; i++) {
		int id = INTEGER(Ractor_ids)[i] - 1;
		if (!(0 <= id && id < size))
			error("actor id out of range");

		dst = matrix_col(&xt, i);
		src = actors_traits(actors, id);

		vector_assign_copy(&dst, src);
	}

	UNPROTECT(1);
	return Rxt;
}

SEXP Riproc_actors_mul(SEXP Ractors, SEXP Rmatrix)
{
	struct actors *actors = Riproc_to_actors(Ractors);
	struct matrix view = Riproc_matrix_view_sexp(Rmatrix);
	int dim = (int)actors_dim(actors);
	int size = (int)actors_size(actors);
	int nrow = (int)matrix_nrow(&view);
	int ncol = (int)matrix_ncol(&view);
	SEXP Rresult;

	if (nrow != dim)
		error("dimension mismatch");

	PROTECT(Rresult = allocMatrix(REALSXP, size, ncol));
	struct matrix result = Riproc_matrix_view_sexp(Rresult);
	struct vector col;
	struct vector dst;

	int j;
	for (j = 0; j < ncol; j++) {
		col = matrix_col(&view, j);
		dst = matrix_col(&result, j);
		actors_mul(1.0, TRANS_NOTRANS, actors, &col, 0.0, &dst);
	}

	UNPROTECT(1);
	return Rresult;
}

SEXP Riproc_actors_tmul(SEXP Ractors, SEXP Rmatrix)
{
	struct actors *actors = Riproc_to_actors(Ractors);
	struct matrix view = Riproc_matrix_view_sexp(Rmatrix);
	int dim = (int)actors_dim(actors);
	int size = (int)actors_size(actors);
	int nrow = (int)matrix_nrow(&view);
	int ncol = (int)matrix_ncol(&view);
	SEXP Rresult;

	if (nrow != size)
		error("dimension mismatch");

	PROTECT(Rresult = allocMatrix(REALSXP, dim, ncol));
	struct matrix result = Riproc_matrix_view_sexp(Rresult);
	struct vector col;
	struct vector dst;

	int j;
	for (j = 0; j < ncol; j++) {
		col = matrix_col(&view, j);
		dst = matrix_col(&result, j);
		actors_mul(1.0, TRANS_TRANS, actors, &col, 0.0, &dst);
	}

	UNPROTECT(1);
	return Rresult;
}
