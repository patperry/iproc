#include "port.h"
#include <assert.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "actors.h"
#include "matrix.h"
#include "r-utils.h"
#include "r-actors.h"

static SEXP Riproc_actors_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_actors_new", (DL_FUNC) & Riproc_actors_new, 1},
	{"Riproc_actors_ngroup", (DL_FUNC) & Riproc_actors_ngroup, 1},
	{"Riproc_actors_size", (DL_FUNC) & Riproc_actors_size, 1},
	{"Riproc_actors_dim", (DL_FUNC) & Riproc_actors_dim, 1},
	{"Riproc_actors_traits", (DL_FUNC) & Riproc_actors_traits, 2},
	{"Riproc_actors_group", (DL_FUNC) & Riproc_actors_group, 2},
	{"Riproc_actors_group_traits", (DL_FUNC) & Riproc_actors_group_traits,
	 2},
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
	iproc_actors *actors = Riproc_to_actors(Ractors);
	iproc_actors_unref(actors);
}

SEXP Riproc_from_actors(iproc_actors * actors)
{
	SEXP Ractors, class;

	iproc_actors_ref(actors);

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

iproc_actors *Riproc_to_actors(SEXP Ractors)
{
	Rboolean null_ok = TRUE;
	SEXP tag = Riproc_actors_type_tag;
	char *type = "actors";
	iproc_actors *actors = Riproc_sexp2ptr(Ractors, null_ok, tag, type);
	return actors;
}

SEXP Riproc_actors_new(SEXP Rtraits_t)
{
	iproc_matrix_view traits_t = Riproc_matrix_view_sexp(Rtraits_t);
	int n = (int)iproc_matrix_ncol(&traits_t.matrix);

	if (!(n > 0))
		error("must specify at least one actor");

	struct vector traits0;
	vector_init_matrix_col(&traits0, &traits_t.matrix, 0);
	ssize_t dim = vector_size(&traits0);

	iproc_actors *actors = iproc_actors_new(dim);

	if (!actors)
		error("could not allocate new actors object");

	struct vector traits;
	int i;
	SEXP Ractors;

	for (i = 0; i < n; i++) {
		vector_init_matrix_col(&traits, &traits_t.matrix, i);
		if (iproc_actors_add(actors, &traits) != i) {
			iproc_actors_unref(actors);
			error
			    ("could not allocate space for %d actors with dim %d",
			     n, dim);
		}
	}

	PROTECT(Ractors = Riproc_from_actors(actors));
	iproc_actors_unref(actors);
	UNPROTECT(1);
	return Ractors;
}

SEXP Riproc_actors_ngroup(SEXP Ractors)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	int ngroup = (int)iproc_actors_ngroup(actors);
	return ScalarInteger(ngroup);
}

SEXP Riproc_actors_size(SEXP Ractors)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	int size = (int)iproc_actors_size(actors);
	return ScalarInteger(size);
}

SEXP Riproc_actors_dim(SEXP Ractors)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	int dim = (int)iproc_actors_dim(actors);
	return ScalarInteger(dim);
}

SEXP Riproc_actors_traits(SEXP Ractors, SEXP Ractor_ids)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	int size = (int)iproc_actors_size(actors);
	int dim = (int)iproc_actors_dim(actors);
	int i, n = GET_LENGTH(Ractor_ids);
	SEXP Rxt;

	PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
	iproc_matrix_view xt = Riproc_matrix_view_sexp(Rxt);
	struct vector dst;
	const struct vector *src;

	for (i = 0; i < n; i++) {
		int id = INTEGER(Ractor_ids)[i] - 1;
		if (!(0 <= id && id < size))
			error("actor id out of range");

		vector_init_matrix_col(&dst, &xt.matrix, i);
		src = iproc_actors_get(actors, id);

		vector_assign_copy(&dst, src);
	}

	UNPROTECT(1);
	return Rxt;
}

SEXP Riproc_actors_group(SEXP Ractors, SEXP Ractor_ids)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	int size = (int)iproc_actors_size(actors);
	int i, n = GET_LENGTH(Ractor_ids);
	SEXP Rgroup_ids;

	PROTECT(Rgroup_ids = NEW_INTEGER(n));

	for (i = 0; i < n; i++) {
		int id = INTEGER(Ractor_ids)[i] - 1;
		if (!(0 <= id && id < size))
			error("actor id out of range");

		int group_id = (int)iproc_actors_group(actors, id);
		INTEGER(Rgroup_ids)[i] = group_id + 1;
	}

	UNPROTECT(1);
	return Rgroup_ids;
}

SEXP Riproc_actors_group_traits(SEXP Ractors, SEXP Rgroup_ids)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	int ngroup = (int)iproc_actors_ngroup(actors);
	int dim = (int)iproc_actors_dim(actors);
	int i, n = GET_LENGTH(Rgroup_ids);
	SEXP Rxt;

	PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
	iproc_matrix_view xt = Riproc_matrix_view_sexp(Rxt);

	struct vector dst;
	const struct vector *src;

	for (i = 0; i < n; i++) {
		int id = INTEGER(Rgroup_ids)[i] - 1;
		if (!(0 <= id && id < ngroup))
			error("group id out of range");

		vector_init_matrix_col(&dst, &xt.matrix, i);
		src = iproc_actors_group_traits(actors, id);

		vector_assign_copy(&dst, src);
	}

	UNPROTECT(1);
	return Rxt;
}

SEXP Riproc_actors_mul(SEXP Ractors, SEXP Rmatrix)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	iproc_matrix_view view = Riproc_matrix_view_sexp(Rmatrix);
	int dim = (int)iproc_actors_dim(actors);
	int size = (int)iproc_actors_size(actors);
	int nrow = (int)iproc_matrix_nrow(&view.matrix);
	int ncol = (int)iproc_matrix_ncol(&view.matrix);
	SEXP Rresult;

	if (nrow != dim)
		error("dimension mismatch");

	PROTECT(Rresult = allocMatrix(REALSXP, size, ncol));
	iproc_matrix_view result = Riproc_matrix_view_sexp(Rresult);
	struct vector col;
	struct vector dst;

	int j;
	for (j = 0; j < ncol; j++) {
		vector_init_matrix_col(&col, &view.matrix, j);
		vector_init_matrix_col(&dst, &result.matrix, j);
		iproc_actors_mul(1.0, IPROC_TRANS_NOTRANS, actors, &col, 0.0,
				 &dst);
	}

	UNPROTECT(1);
	return Rresult;
}

SEXP Riproc_actors_tmul(SEXP Ractors, SEXP Rmatrix)
{
	iproc_actors *actors = Riproc_to_actors(Ractors);
	iproc_matrix_view view = Riproc_matrix_view_sexp(Rmatrix);
	int dim = (int)iproc_actors_dim(actors);
	int size = (int)iproc_actors_size(actors);
	int nrow = (int)iproc_matrix_nrow(&view.matrix);
	int ncol = (int)iproc_matrix_ncol(&view.matrix);
	SEXP Rresult;

	if (nrow != size)
		error("dimension mismatch");

	PROTECT(Rresult = allocMatrix(REALSXP, dim, ncol));
	iproc_matrix_view result = Riproc_matrix_view_sexp(Rresult);
	struct vector col;
	struct vector dst;

	int j;
	for (j = 0; j < ncol; j++) {
		vector_init_matrix_col(&col, &view.matrix, j);
		vector_init_matrix_col(&dst, &result.matrix, j);
		iproc_actors_mul(1.0, IPROC_TRANS_TRANS, actors, &col, 0.0,
				 &dst);
	}

	UNPROTECT(1);
	return Rresult;
}
