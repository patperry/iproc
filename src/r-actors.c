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
    { "Riproc_actors_vector", (DL_FUNC) &Riproc_actors_vector, 2 },
    { "Riproc_actors_class",  (DL_FUNC) &Riproc_actors_class,  2 },
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

SEXP
Riproc_actors_vector (SEXP Ractors,
                      SEXP Ractor_ids)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t size = iproc_actors_size(actors);
    int64_t dim = iproc_actors_dim(actors);
    int64_t i, n = GET_LENGTH(Ractor_ids);
    SEXP *Rxt;

    PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
    double *data = NUMERIC_POINTER(Rxt);
    iproc_matrix_view xt = iproc_matrix_view_array(data, dim, n);

    for (i = 0; i < n; i++) {
        int64_t id = INTEGER(Ractor_ids)[i] - 1;
        if (!(0 <= id && id < size))
            error("actor id out of range");

        iproc_vector_view dst = iproc_matrix_col(&xt.matrix, i);
        iproc_vector *src = iproc_actors_vector(actors, id);

        iproc_vector_copy(&dst.vector, src);
    }

    UNPROTECT(1);
    return Rxt;
}


SEXP
Riproc_actors_class (SEXP Ractors,
                     SEXP Ractor_ids)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t size = iproc_actors_size(actors);
    int64_t i, n = GET_LENGTH(Ractor_ids);
    SEXP *Rclass_ids;

    PROTECT(Rclass_ids = NEW_INTEGER(n));

    for (i = 0; i < n; i++) {
        int64_t id = INTEGER(Ractor_ids)[i] - 1;
        if (!(0 <= id && id < size))
            error("actor id out of range");

        int64_t class_id = iproc_actors_class(actors, id);
        INTEGER(Rclass_ids)[i] = class_id + 1;
    }

    UNPROTECT(1);
    return Rclass_ids;
}


SEXP
Riproc_actors_class_vector (SEXP Ractors,
                            SEXP Rclass_ids)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t nclass = iproc_actors_nclass(actors);
    int64_t dim = iproc_actors_dim(actors);
    int64_t i, n = GET_LENGTH(Rclass_ids);
    SEXP *Rxt;

    PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
    double *data = NUMERIC_POINTER(Rxt);
    iproc_matrix_view xt = iproc_matrix_view_array(data, dim, n);

    for (i = 0; i < n; i++) {
        int64_t id = INTEGER(Rclass_ids)[i] - 1;
        if (!(0 <= id && id < nclass))
            error("class id out of range");

        iproc_vector_view dst = iproc_matrix_col(&xt.matrix, i);
        iproc_vector *src = iproc_actors_class_vector(actors, id);

        iproc_vector_copy(&dst.vector, src);
    }

    UNPROTECT(1);
    return Rxt;
}
