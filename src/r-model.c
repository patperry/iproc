#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-actors.h"
#include "r-vars.h"
#include "r-model.h"


static SEXP Riproc_model_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_model_new",       (DL_FUNC) &Riproc_model_new,       3 },
    { "Riproc_model_vars",      (DL_FUNC) &Riproc_model_vars,      1 },
    { "Riproc_model_coefs",     (DL_FUNC) &Riproc_model_coefs,     1 },
    { "Riproc_model_has_loops", (DL_FUNC) &Riproc_model_has_loops, 1 },
    { "Riproc_model_dim",       (DL_FUNC) &Riproc_model_dim,       1 },
    { "Riproc_model_nreceiver", (DL_FUNC) &Riproc_model_nreceiver, 1 },
    { "Riproc_model_nsender",   (DL_FUNC) &Riproc_model_nsender,   1 },
    { NULL,                     NULL,                              0 }
};


void
Riproc_model_init (DllInfo *info)
{
    Riproc_model_type_tag = install("Riproc_model_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_model_free (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    iproc_model_unref(model);
}


iproc_model *
Riproc_to_model (SEXP Rmodel)
{
    iproc_model *model = Riproc_sexp2ptr(Rmodel, FALSE, Riproc_model_type_tag, "model");
    return model;
}

SEXP
Riproc_from_model (iproc_model *model)
{
    SEXP Rmodel, class;

    iproc_model_ref(model);

    Rmodel = R_MakeExternalPtr(model, Riproc_model_type_tag, R_NilValue);
    R_RegisterCFinalizer(Rmodel, Riproc_model_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("model"));
    classgets(Rmodel, class);
    UNPROTECT(1);

    return Rmodel;
}

SEXP
Riproc_model_new (SEXP Rvars,
                  SEXP Rcoefs,
                  SEXP Rhas_loops)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    iproc_vector_view coefs = Riproc_vector_view_sexp(Rcoefs);
    Rboolean has_loops = LOGICAL_VALUE(Rhas_loops);

    if (iproc_vars_dim(vars) != iproc_vector_dim(&coefs.vector))
        error("vars and coefs have different dimensions");
    if (has_loops == NA_LOGICAL)
        error("has.loops is be NA");

    iproc_model *model = iproc_model_new(vars, &coefs.vector, has_loops);
    SEXP Rmodel = Riproc_from_model(model);
    iproc_model_unref(model);
    return Rmodel;
}

SEXP
Riproc_model_dim (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    int64_t dim = iproc_model_dim(model);
    return ScalarInteger(dim);
}

SEXP
Riproc_model_nsender (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    int64_t n = iproc_model_nsender(model);
    return ScalarInteger(n);
}

SEXP
Riproc_model_nreceiver (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    int64_t n = iproc_model_nreceiver(model);
    return ScalarInteger(n);
}

SEXP
Riproc_model_vars (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    iproc_vars *vars = iproc_model_vars(model);
    return Riproc_from_vars(vars);
}

SEXP
Riproc_model_coefs (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    iproc_vector *coefs = iproc_model_coefs(model);

    return Riproc_vector_new_copy(coefs);
}

SEXP
Riproc_model_has_loops (SEXP Rmodel)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    if (iproc_model_has_loops(model)) {
        return ScalarLogical(TRUE);
    } else {
        return ScalarLogical(FALSE);
    }
}
