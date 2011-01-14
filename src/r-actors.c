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
    { "Riproc_actors_new",           (DL_FUNC) &Riproc_actors_new,           2 },
    { "Riproc_actors_ngroup",        (DL_FUNC) &Riproc_actors_ngroup,        1 },
    { "Riproc_actors_size",          (DL_FUNC) &Riproc_actors_size,          1 },
    { "Riproc_actors_dim",           (DL_FUNC) &Riproc_actors_dim,           1 },
    { "Riproc_actors_traits",        (DL_FUNC) &Riproc_actors_traits,        2 },
    { "Riproc_actors_group",         (DL_FUNC) &Riproc_actors_group,         2 },
    { "Riproc_actors_group_traits",  (DL_FUNC) &Riproc_actors_group_traits,  2 },
    { NULL,                          NULL,                                   0 }
};


void
Riproc_actors_init (DllInfo *info)
{
    Riproc_actors_type_tag = install("Riproc_actors_type_tag");
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

    iproc_actors_ref(actors);

    /* store the actors pointer in an external pointer */
    PROTECT(Ractors = R_MakeExternalPtr(actors, Riproc_actors_type_tag, R_NilValue));
    R_RegisterCFinalizer(Ractors, Riproc_actors_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("actors"));
    classgets(Ractors, class);

    UNPROTECT(2);
    return Ractors;
}

iproc_actors *
Riproc_to_actors (SEXP Ractors)
{
    Rboolean null_ok = TRUE;
    SEXP tag = Riproc_actors_type_tag;
    char *type = "actors";
    iproc_actors *actors = Riproc_sexp2ptr(Ractors, null_ok, tag, type);
    return actors;
}

SEXP
Riproc_actors_new (SEXP Rgroups,
                   SEXP Rgroup_traits_t)
{
    int64_t n = GET_LENGTH(Rgroups);
    int *groups = INTEGER_POINTER(Rgroups);
    iproc_matrix_view group_traits_t = Riproc_matrix_view_sexp(Rgroup_traits_t);
    int64_t k = iproc_matrix_ncol(&group_traits_t.matrix);

    if (!(k > 0))
        error("must specify at least one group");

    iproc_vector_view traits0 = iproc_matrix_col(&group_traits_t.matrix, 0);
    iproc_actors *actors = iproc_actors_new(n, &traits0.vector);
    int64_t i;
    SEXP Ractors;

    for (i = 1; i < k; i++) {
        iproc_vector_view traits = iproc_matrix_col(&group_traits_t.matrix, i);
        iproc_actors_append_group(actors, &traits.vector);
    }

    for (i = 0; i < n; i++) {
        int64_t c = groups[i] - 1;
        
        if (!(0 <= c && c < k))
            error("group id out of bounds");

        iproc_actors_set(actors, i, c);
    }

    PROTECT(Ractors = Riproc_from_actors(actors));
    iproc_actors_unref(actors);

    UNPROTECT(1);
    return Ractors;
}

SEXP
Riproc_actors_ngroup (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t ngroup = iproc_actors_ngroup(actors);
    return ScalarInteger(ngroup);
}

SEXP
Riproc_actors_size (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t size = iproc_actors_size(actors);
    return ScalarInteger(size);
}

SEXP
Riproc_actors_dim (SEXP Ractors)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t dim = iproc_actors_dim(actors);
    return ScalarInteger(dim);
}

SEXP
Riproc_actors_traits (SEXP Ractors,
                      SEXP Ractor_ids)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t size = iproc_actors_size(actors);
    int64_t dim = iproc_actors_dim(actors);
    int64_t i, n = GET_LENGTH(Ractor_ids);
    SEXP Rxt;

    PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
    iproc_matrix_view xt = Riproc_matrix_view_sexp(Rxt);

    for (i = 0; i < n; i++) {
        int64_t id = INTEGER(Ractor_ids)[i] - 1;
        if (!(0 <= id && id < size))
            error("actor id out of range");

        iproc_vector_view dst = iproc_matrix_col(&xt.matrix, i);
        iproc_vector *src = iproc_actors_traits(actors, id);

        iproc_vector_copy(&dst.vector, src);
    }

    UNPROTECT(1);
    return Rxt;
}


SEXP
Riproc_actors_group (SEXP Ractors,
                     SEXP Ractor_ids)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t size = iproc_actors_size(actors);
    int64_t i, n = GET_LENGTH(Ractor_ids);
    SEXP Rgroup_ids;

    PROTECT(Rgroup_ids = NEW_INTEGER(n));

    for (i = 0; i < n; i++) {
        int64_t id = INTEGER(Ractor_ids)[i] - 1;
        if (!(0 <= id && id < size))
            error("actor id out of range");

        int64_t group_id = iproc_actors_group(actors, id);
        INTEGER(Rgroup_ids)[i] = group_id + 1;
    }

    UNPROTECT(1);
    return Rgroup_ids;
}


SEXP
Riproc_actors_group_traits (SEXP Ractors,
                            SEXP Rgroup_ids)
{
    iproc_actors *actors = Riproc_to_actors(Ractors);
    int64_t ngroup = iproc_actors_ngroup(actors);
    int64_t dim = iproc_actors_dim(actors);
    int64_t i, n = GET_LENGTH(Rgroup_ids);
    SEXP Rxt;

    PROTECT(Rxt = allocMatrix(REALSXP, dim, n));
    iproc_matrix_view xt = Riproc_matrix_view_sexp(Rxt);

    for (i = 0; i < n; i++) {
        int64_t id = INTEGER(Rgroup_ids)[i] - 1;
        if (!(0 <= id && id < ngroup))
            error("group id out of range");

        iproc_vector_view dst = iproc_matrix_col(&xt.matrix, i);
        iproc_vector *src = iproc_actors_group_traits(actors, id);

        iproc_vector_copy(&dst.vector, src);
    }

    UNPROTECT(1);
    return Rxt;
}
