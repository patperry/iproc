#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-actors.h"
#include "r-cursor.h"
#include "r-design.h"


static SEXP Riproc_design_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_design_new",       (DL_FUNC) &Riproc_design_new,       3 },
    { "Riproc_design_dim",       (DL_FUNC) &Riproc_design_dim,       1 },
    { "Riproc_design_nreceiver", (DL_FUNC) &Riproc_design_nreceiver, 1 },
    { "Riproc_design_nsender",   (DL_FUNC) &Riproc_design_nsender,   1 },
    { "Riproc_design_receivers", (DL_FUNC) &Riproc_design_receivers, 1 },
    { "Riproc_design_senders",   (DL_FUNC) &Riproc_design_senders,   1 },
    { "Riproc_design_mul",       (DL_FUNC) &Riproc_design_mul,       4 },
    { "Riproc_design_tmul",      (DL_FUNC) &Riproc_design_tmul,      4 },
    { NULL,                    NULL,                             0 }
};


void
Riproc_design_init (DllInfo *info)
{
    Riproc_design_type_tag = install("Riproc_design_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_design_free (SEXP Rdesign)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    iproc_design_unref(design);
}


iproc_design *
Riproc_to_design (SEXP Rdesign)
{
    iproc_design *design =  Riproc_sexp2ptr(Rdesign, TRUE, Riproc_design_type_tag, "design");
    return design;
}

SEXP
Riproc_from_design (iproc_design *design)
{
    SEXP Rdesign, class;

    iproc_design_ref(design);

    PROTECT(Rdesign = R_MakeExternalPtr(design, Riproc_design_type_tag, R_NilValue));
    R_RegisterCFinalizer(Rdesign, Riproc_design_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("iproc.design"));
    classgets(Rdesign, class);

    UNPROTECT(2);
    return Rdesign;
}

SEXP
Riproc_design_new (SEXP Rsenders,
                 SEXP Rreceivers,
                 SEXP Rrecip_intervals)
{
    iproc_actors *senders = Riproc_to_actors(Rsenders);
    iproc_actors *receivers = Riproc_to_actors(Rreceivers);
    Riproc_design_udata *udata = Riproc_design_udata_new(Rrecip_intervals);
    int64_t ndynamic = Riproc_design_udata_dim(udata);
    iproc_design *design = iproc_design_new(senders, receivers, ndynamic, udata,
                                      Riproc_design_udata_get_sdesign_vars,
                                      Riproc_design_udata_free);
    SEXP Rdesign;

    PROTECT(Rdesign = Riproc_from_design(design));
    iproc_design_unref(design);

    UNPROTECT(1);
    return Rdesign;
}

SEXP
Riproc_design_dim (SEXP Rdesign)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    int dim = (int)iproc_design_dim(design);
    return ScalarInteger(dim);
}


SEXP
Riproc_design_nsender (SEXP Rdesign)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    int nsender = (int)iproc_design_nsender(design);
    return ScalarInteger(nsender);
}

SEXP
Riproc_design_nreceiver (SEXP Rdesign)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    int nreceiver = (int)iproc_design_nreceiver(design);
    return ScalarInteger(nreceiver);
}

SEXP
Riproc_design_senders (SEXP Rdesign)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    iproc_actors *senders = iproc_design_senders(design);
    return Riproc_from_actors(senders);
}

SEXP
Riproc_design_receivers (SEXP Rdesign)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    iproc_actors *receivers = iproc_design_receivers(design);
    return Riproc_from_actors(receivers);
}

SEXP
Riproc_design_mul (SEXP Rdesign,
                 SEXP Rx,
                 SEXP Rsender,
                 SEXP Rcursor)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    int dim = (int)iproc_design_dim(design);
    int nsender = (int)iproc_design_nsender(design);
    int nreceiver = (int)iproc_design_nreceiver(design);
    iproc_matrix_view x = Riproc_matrix_view_sexp(Rx);
    int nrow = (int)iproc_matrix_nrow(&x.matrix);
    int ncol = (int)iproc_matrix_ncol(&x.matrix);
    int sender = INTEGER(Rsender)[0] - 1;
    iproc_message_iter *cursor =  (Rcursor == NULL_USER_OBJECT
                                   ? NULL
                                   : Riproc_to_cursor(Rcursor));
    iproc_history *history = iproc_message_iter_history(cursor);
    
    if (sender < 0 || sender >= nsender)
        error("invalid sender");
    if (nrow != dim)
        error("dimension mismatch");

    SEXP Rresult;
    PROTECT(Rresult = allocMatrix(REALSXP, nreceiver, ncol));
    iproc_matrix_view result = Riproc_matrix_view_sexp(Rresult);
    iproc_design_ctx *ctx = iproc_design_ctx_new(design, sender, history);

    int j;
    for (j = 0; j < ncol; j++) {
        iproc_vector_view col = iproc_matrix_col(&x.matrix, j);
        iproc_vector_view dst = iproc_matrix_col(&result.matrix, j);
        iproc_design_ctx_mul(1.0, IPROC_TRANS_NOTRANS, ctx, &col.vector, 0.0, &dst.vector);
    }

    iproc_design_ctx_unref(ctx);

    UNPROTECT(1);
    return Rresult;
}

SEXP
Riproc_design_tmul (SEXP Rdesign,
                  SEXP Rx,
                  SEXP Rsender,
                  SEXP Rcursor)
{
    iproc_design *design = Riproc_to_design(Rdesign);
    int dim = (int)iproc_design_dim(design);
    int nsender = (int)iproc_design_nsender(design);
    int nreceiver = (int)iproc_design_nreceiver(design);
    iproc_matrix_view x = Riproc_matrix_view_sexp(Rx);
    int nrow = (int)iproc_matrix_nrow(&x.matrix);
    int ncol = (int)iproc_matrix_ncol(&x.matrix);
    int sender = INTEGER(Rsender)[0] - 1;
    iproc_message_iter *cursor =  (Rcursor == NULL_USER_OBJECT
                                   ? NULL
                                   : Riproc_to_cursor(Rcursor));
    iproc_history *history = iproc_message_iter_history(cursor);
    
    if (sender < 0 || sender >= nsender)
        error("invalid sender");
    if (nrow != nreceiver)
        error("dimension mismatch");

    SEXP Rresult;
    PROTECT(Rresult = allocMatrix(REALSXP, dim, ncol));
    iproc_matrix_view result = Riproc_matrix_view_sexp(Rresult);
    iproc_design_ctx *ctx = iproc_design_ctx_new(design, sender, history);

    int j;
    for (j = 0; j < ncol; j++) {
        iproc_vector_view col = iproc_matrix_col(&x.matrix, j);
        iproc_vector_view dst = iproc_matrix_col(&result.matrix, j);
        iproc_design_ctx_mul(1.0, IPROC_TRANS_TRANS, ctx, &col.vector, 0.0, &dst.vector);
    }

    iproc_design_ctx_unref(ctx);

    UNPROTECT(1);
    return Rresult;
}