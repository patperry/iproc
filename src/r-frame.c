#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-actors.h"
#include "r-cursor.h"
#include "r-frame.h"


static SEXP Riproc_frame_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_frame_new",       (DL_FUNC) &Riproc_frame_new,       3 },
    { "Riproc_frame_dim",       (DL_FUNC) &Riproc_frame_dim,       1 },
    { "Riproc_frame_nreceiver", (DL_FUNC) &Riproc_frame_nreceiver, 1 },
    { "Riproc_frame_nsender",   (DL_FUNC) &Riproc_frame_nsender,   1 },
    { "Riproc_frame_receivers", (DL_FUNC) &Riproc_frame_receivers, 1 },
    { "Riproc_frame_senders",   (DL_FUNC) &Riproc_frame_senders,   1 },
    { "Riproc_frame_mul",       (DL_FUNC) &Riproc_frame_mul,       4 },
    { "Riproc_frame_tmul",      (DL_FUNC) &Riproc_frame_tmul,      4 },
    { NULL,                    NULL,                             0 }
};


void
Riproc_frame_init (DllInfo *info)
{
    Riproc_frame_type_tag = install("Riproc_frame_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_frame_free (SEXP Rframe)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    iproc_frame_unref(frame);
}


iproc_frame *
Riproc_to_frame (SEXP Rframe)
{
    iproc_frame *frame =  Riproc_sexp2ptr(Rframe, TRUE, Riproc_frame_type_tag, "frame");
    return frame;
}

SEXP
Riproc_from_frame (iproc_frame *frame)
{
    SEXP Rframe, class;

    iproc_frame_ref(frame);

    PROTECT(Rframe = R_MakeExternalPtr(frame, Riproc_frame_type_tag, R_NilValue));
    R_RegisterCFinalizer(Rframe, Riproc_frame_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("iproc.frame"));
    classgets(Rframe, class);

    UNPROTECT(2);
    return Rframe;
}

SEXP
Riproc_frame_new (SEXP Rsenders,
                 SEXP Rreceivers,
                 SEXP Rreceive_intervals)
{
    iproc_actors *senders = Riproc_to_actors(Rsenders);
    iproc_actors *receivers = Riproc_to_actors(Rreceivers);
    Riproc_frame_udata *udata = Riproc_frame_udata_new(Rreceive_intervals);
    int64_t ndynamic = Riproc_frame_udata_dim(udata);
    iproc_frame *frame = iproc_frame_new(senders, receivers, ndynamic, udata,
                                      Riproc_frame_udata_get_sender_frame,
                                      Riproc_frame_udata_free);
    SEXP Rframe;

    PROTECT(Rframe = Riproc_from_frame(frame));
    iproc_frame_unref(frame);

    UNPROTECT(1);
    return Rframe;
}

SEXP
Riproc_frame_dim (SEXP Rframe)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    int64_t dim = iproc_frame_dim(frame);
    return ScalarInteger(dim);
}


SEXP
Riproc_frame_nsender (SEXP Rframe)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    int64_t nsender = iproc_frame_nsender(frame);
    return ScalarInteger(nsender);
}

SEXP
Riproc_frame_nreceiver (SEXP Rframe)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    int64_t nreceiver = iproc_frame_nreceiver(frame);
    return ScalarInteger(nreceiver);
}

SEXP
Riproc_frame_senders (SEXP Rframe)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    iproc_actors *senders = iproc_frame_senders(frame);
    return Riproc_from_actors(senders);
}

SEXP
Riproc_frame_receivers (SEXP Rframe)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    iproc_actors *receivers = iproc_frame_receivers(frame);
    return Riproc_from_actors(receivers);
}

SEXP
Riproc_frame_mul (SEXP Rframe,
                 SEXP Rx,
                 SEXP Rsender,
                 SEXP Rcursor)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    int64_t dim = iproc_frame_dim(frame);
    int64_t nsender = iproc_frame_nsender(frame);
    int64_t nreceiver = iproc_frame_nreceiver(frame);
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
    iproc_frame_ctx *ctx = iproc_frame_ctx_new(frame, sender, history);

    int64_t j;
    for (j = 0; j < ncol; j++) {
        iproc_vector_view col = iproc_matrix_col(&x.matrix, j);
        iproc_vector_view dst = iproc_matrix_col(&result.matrix, j);
        iproc_frame_ctx_mul(1.0, IPROC_TRANS_NOTRANS, ctx, &col.vector, 0.0, &dst.vector);
    }

    iproc_frame_ctx_unref(ctx);

    UNPROTECT(1);
    return Rresult;
}

SEXP
Riproc_frame_tmul (SEXP Rframe,
                  SEXP Rx,
                  SEXP Rsender,
                  SEXP Rcursor)
{
    iproc_frame *frame = Riproc_to_frame(Rframe);
    int64_t dim = iproc_frame_dim(frame);
    int64_t nsender = iproc_frame_nsender(frame);
    int64_t nreceiver = iproc_frame_nreceiver(frame);
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
    iproc_frame_ctx *ctx = iproc_frame_ctx_new(frame, sender, history);

    int64_t j;
    for (j = 0; j < ncol; j++) {
        iproc_vector_view col = iproc_matrix_col(&x.matrix, j);
        iproc_vector_view dst = iproc_matrix_col(&result.matrix, j);
        iproc_frame_ctx_mul(1.0, IPROC_TRANS_TRANS, ctx, &col.vector, 0.0, &dst.vector);
    }

    iproc_frame_ctx_unref(ctx);

    UNPROTECT(1);
    return Rresult;
}
