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
    { "Riproc_actors_new",    (DL_FUNC) &Riproc_actors_new,    2 },
    { "Riproc_actors_nclass", (DL_FUNC) &Riproc_actors_nclass, 1 },
    { "Riproc_actors_size",   (DL_FUNC) &Riproc_actors_size,   1 },
    { "Riproc_actors_dim",    (DL_FUNC) &Riproc_actors_dim,    1 },
    { NULL,                   NULL,                            0 }
};


void
Riproc_actors_init (DllInfo *info)
{
    Riproc_actors_type_tag = install("Riproc_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_actors_free (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    iproc_actors_unref(actors);
}

SEXP
Riproc_from_actors (iproc_actors *actors)
{
    SEXP Ractors, class;
    iproc_actors *actors1;

    /* store the actors pointer in an external pointer */
    Ractors = R_MakeExternalPtr(actors, Riproc_actors_type_tag, R_NilValue);
    R_RegisterCFinalizer(Ractors, Riproc_actors_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("actors"));
    classgets(Ractors, class);
    UNPROTECT(1);

    return Ractors;
}

iproc_actors *
Riproc_to_actors (SEXP Ractors)
{
    Rboolean null_ok = TRUE;
    SEXP tag = Riproc_actors_type_tag;
    char *type = "actors";

    return Riproc_sexp2ptr(Ractors, null_ok, tag, type);
}

SEXP
Riproc_actors_new (SEXP Rclasses,
                   SEXP Rclass_vectors_t)
{
    int64_t n = GET_LENGTH(Rclasses);
    int *classes = INTEGER_POINTER(Rclasses);
    int64_t p    = INTEGER(GET_DIM(Rclass_vectors_t))[0];
    int64_t k    = INTEGER(GET_DIM(Rclass_vectors_t))[1];
    double *data = NUMERIC_POINTER(Rclass_vectors_t);
    iproc_matrix_view class_vectors_t = iproc_matrix_view_array(data, p, k);
    iproc_vector_view defvector = iproc_matrix_col(&class_vectors_t.matrix, 0);
    iproc_actors *actors = iproc_actors_new(n, &defvector.vector);
    int64_t i;
    SEXP Ractors;

    for (i = 1; i < k; i++) {
        iproc_vector_view col = iproc_matrix_col(&class_vectors_t.matrix, i);
        iproc_actors_append_class(actors, &col.vector);
    }

    for (i = 0; i < n; i++) {
        int64_t c = classes[i] - 1;
        iproc_actors_set(actors, i, c);
    }

    printf("Allocated actors %p\n", actors);
    Ractors = Riproc_from_actors(actors);
    return Ractors;
}

SEXP
Riproc_actors_nclass (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t nclass = iproc_actors_nclass(actors);
    SEXP Rnclass;

    PROTECT(Rnclass = NEW_INTEGER(1));
    INTEGER(Rnclass)[0] = nclass;

    UNPROTECT(1);
    return Rnclass;
}

SEXP
Riproc_actors_size (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t size = iproc_actors_size(actors);
    SEXP Rsize;

    PROTECT(Rsize = NEW_INTEGER(1));
    INTEGER(Rsize)[0] = size;

    UNPROTECT(1);
    return Rsize;
}

SEXP
Riproc_actors_dim (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t dim = iproc_actors_dim(actors);
    SEXP Rdim;

    PROTECT(Rdim = NEW_INTEGER(1));
    INTEGER(Rdim)[0] = dim;

    UNPROTECT(1);
    return Rdim;
}
