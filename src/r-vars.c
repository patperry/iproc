#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-actors.h"
#include "r-cursor.h"
#include "r-vars.h"


static SEXP Riproc_vars_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_vars_new",       (DL_FUNC) &Riproc_vars_new,       2 },
    { "Riproc_vars_dim",       (DL_FUNC) &Riproc_vars_dim,       1 },
    { "Riproc_vars_nreceiver", (DL_FUNC) &Riproc_vars_nreceiver, 1 },
    { "Riproc_vars_nsender",   (DL_FUNC) &Riproc_vars_nsender,   1 },
    { "Riproc_vars_receivers", (DL_FUNC) &Riproc_vars_receivers, 1 },
    { "Riproc_vars_senders",   (DL_FUNC) &Riproc_vars_senders,   1 },
    { "Riproc_vars_mul",       (DL_FUNC) &Riproc_vars_mul,       4 },
    { "Riproc_vars_tmul",      (DL_FUNC) &Riproc_vars_tmul,      4 },
    { NULL,                    NULL,                             0 }
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

    PROTECT(Rvars = R_MakeExternalPtr(vars, Riproc_vars_type_tag, R_NilValue));
    R_RegisterCFinalizer(Rvars, Riproc_vars_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("vars"));
    classgets(Rvars, class);

    UNPROTECT(2);
    return Rvars;
}

SEXP
Riproc_vars_new (SEXP Rsenders,
                 SEXP Rreceivers)
{
    iproc_actors *senders = Riproc_to_actors(Rsenders);
    iproc_actors *receivers = Riproc_to_actors(Rreceivers);
    int64_t ndynamic = 0;
    iproc_vars *vars = iproc_vars_new(senders, receivers, ndynamic, NULL);
    SEXP Rvars;

    PROTECT(Rvars = Riproc_from_vars(vars));
    iproc_vars_unref(vars);

    UNPROTECT(1);
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

SEXP
Riproc_vars_mul (SEXP Rvars,
                 SEXP Rx,
                 SEXP Rsender,
                 SEXP Rcursor)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    int64_t dim = iproc_vars_dim(vars);
    int64_t nsender = iproc_vars_nsender(vars);
    int64_t nreceiver = iproc_vars_nreceiver(vars);
    iproc_matrix_view x = Riproc_matrix_view_sexp(Rx);
    int64_t nrow = iproc_matrix_nrow(&x.matrix);
    int64_t ncol = iproc_matrix_ncol(&x.matrix);
    int64_t sender = INTEGER(Rsender)[0] - 1;
    iproc_cursor *cursor =  (Rcursor == NULL_USER_OBJECT
                            ? NULL
                            : Riproc_to_cursor(Rcursor));
    iproc_history *history = iproc_cursor_history(cursor);
    
    if (sender < 0 || sender >= nsender)
        error("invalid sender");
    if (nrow != dim)
        error("dimension mismatch");

    SEXP Rresult;
    PROTECT(Rresult = allocMatrix(REALSXP, nreceiver, ncol));
    iproc_matrix_view result = Riproc_matrix_view_sexp(Rresult);
    iproc_vars_ctx *ctx = iproc_vars_ctx_new(vars, sender, history);

    int64_t j;
    for (j = 0; j < ncol; j++) {
        iproc_vector_view col = iproc_matrix_col(&x.matrix, j);
        iproc_vector_view dst = iproc_matrix_col(&result.matrix, j);
        iproc_vars_ctx_mul(1.0, IPROC_TRANS_NOTRANS, ctx, &col.vector, 0.0, &dst.vector);
    }

    iproc_vars_ctx_unref(ctx);

    UNPROTECT(1);
    return Rresult;
}

SEXP
Riproc_vars_tmul (SEXP Rvars,
                  SEXP Rx,
                  SEXP Rsender,
                  SEXP Rcursor)
{
    iproc_vars *vars = Riproc_to_vars(Rvars);
    int64_t dim = iproc_vars_dim(vars);
    int64_t nsender = iproc_vars_nsender(vars);
    int64_t nreceiver = iproc_vars_nreceiver(vars);
    iproc_matrix_view x = Riproc_matrix_view_sexp(Rx);
    int64_t nrow = iproc_matrix_nrow(&x.matrix);
    int64_t ncol = iproc_matrix_ncol(&x.matrix);
    int64_t sender = INTEGER(Rsender)[0] - 1;
    iproc_cursor *cursor =  (Rcursor == NULL_USER_OBJECT
                            ? NULL
                            : Riproc_to_cursor(Rcursor));
    iproc_history *history = iproc_cursor_history(cursor);
    
    if (sender < 0 || sender >= nsender)
        error("invalid sender");
    if (nrow != nreceiver)
        error("dimension mismatch");

    SEXP Rresult;
    PROTECT(Rresult = allocMatrix(REALSXP, dim, ncol));
    iproc_matrix_view result = Riproc_matrix_view_sexp(Rresult);
    iproc_vars_ctx *ctx = iproc_vars_ctx_new(vars, sender, history);

    int64_t j;
    for (j = 0; j < ncol; j++) {
        iproc_vector_view col = iproc_matrix_col(&x.matrix, j);
        iproc_vector_view dst = iproc_matrix_col(&result.matrix, j);
        iproc_vars_ctx_mul(1.0, IPROC_TRANS_TRANS, ctx, &col.vector, 0.0, &dst.vector);
    }

    iproc_vars_ctx_unref(ctx);

    UNPROTECT(1);
    return Rresult;
}
