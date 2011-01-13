#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-actors.h"
#include "r-vars.h"


static SEXP Riproc_vars_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_vars_new",       (DL_FUNC) &Riproc_vars_new,       2 },
    { "Riproc_vars_dim",       (DL_FUNC) &Riproc_vars_dim,       1 },
    { "Riproc_vars_nreceiver", (DL_FUNC) &Riproc_vars_nreceiver, 1 },
    { "Riproc_vars_nsender",   (DL_FUNC) &Riproc_vars_nsender,   1 },
    { "Riproc_vars_receivers", (DL_FUNC) &Riproc_vars_receivers, 1 },
    { "Riproc_vars_senders",   (DL_FUNC) &Riproc_vars_senders,   1 },
    { NULL,                          NULL,                       0 }
};


void
Riproc_vars_init (DllInfo *info)
{
    Riproc_vars_type_tag = install("Riproc_vars_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_vars_free (SEXP Rvars)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    iproc_vars_unref(vars);
}


iproc_vars *
Riproc_to_vars (SEXP Rvars)
{
    iproc_vars *vars =  Riproc_sexp2ptr(Rvars, TRUE, Riproc_vars_type_tag, "vars");
    return vars;
}

SEXP
Riproc_from_vars (iproc_vars *vars)
{
    SEXP Rvars, class;

    iproc_vars_ref(vars);

    Rvars = R_MakeExternalPtr(vars, Riproc_vars_type_tag, R_NilValue);
    R_RegisterCFinalizer(Rvars, Riproc_vars_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("vars"));
    classgets(Rvars, class);
    UNPROTECT(1);

    return Rvars;
}

SEXP
Riproc_vars_new (SEXP Rsenders,
                 SEXP Rreceivers)
{
    iproc_actors *senders = Riproc_to_actors(Rsenders);
    iproc_actors *receivers = Riproc_to_actors(Rreceivers);
    iproc_vars *vars = iproc_vars_new(senders, receivers);
    SEXP Rvars = Riproc_from_vars(vars);
    iproc_vars_unref(vars);
    return Rvars;
}

SEXP
Riproc_vars_dim (SEXP Rvars)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    int64_t dim = iproc_vars_dim(vars);
    return ScalarInteger(dim);
}


SEXP
Riproc_vars_nsender (SEXP Rvars)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    int64_t nsender = iproc_vars_nsender(vars);
    return ScalarInteger(nsender);
}

SEXP
Riproc_vars_nreceiver (SEXP Rvars)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    int64_t nreceiver = iproc_vars_nreceiver(vars);
    return ScalarInteger(nreceiver);
}

SEXP
Riproc_vars_senders (SEXP Rvars)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    iproc_actors *senders = iproc_vars_senders(vars);
    return Riproc_from_actors(senders);
}

SEXP
Riproc_vars_receivers (SEXP Rvars)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    iproc_actors *receivers = iproc_vars_receivers(vars);
    return Riproc_from_actors(receivers);
}
