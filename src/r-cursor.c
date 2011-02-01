

#include "r-messages.h"
#include "r-utils.h"
#include "r-cursor.h"



static SEXP Riproc_cursor_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_cursor_new",      (DL_FUNC) &Riproc_cursor_new,      2 },
    { "Riproc_cursor_advance",  (DL_FUNC) &Riproc_cursor_advance,  1 },
    { "Riproc_cursor_started",  (DL_FUNC) &Riproc_cursor_started,  1 },
    { "Riproc_cursor_finished", (DL_FUNC) &Riproc_cursor_finished, 1 },
    { "Riproc_cursor_reset",    (DL_FUNC) &Riproc_cursor_reset,    1 },
    { "Riproc_cursor_time",     (DL_FUNC) &Riproc_cursor_time,     1 },
    { "Riproc_cursor_nties",    (DL_FUNC) &Riproc_cursor_nties,    1 },
    { "Riproc_cursor_from",     (DL_FUNC) &Riproc_cursor_from,     1 },
    { "Riproc_cursor_to",       (DL_FUNC) &Riproc_cursor_to,       1 },
    { NULL,                     NULL,                               0 }
};

void
Riproc_cursor_init (DllInfo *info)
{
    Riproc_cursor_type_tag = install("Riproc_cursor_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_cursor_free (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);
    iproc_cursor_unref(cursor);
}


SEXP
Riproc_from_cursor (iproc_cursor *cursor)
{
    SEXP Rcursor, class;

    iproc_cursor_ref(cursor);

    /* store the cursor pointer in an external pointer */
    PROTECT(Rcursor = R_MakeExternalPtr(cursor, Riproc_cursor_type_tag, R_NilValue));
    R_RegisterCFinalizer(Rcursor, Riproc_cursor_free);

    /* set the class of the result */
    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("cursor"));
    classgets(Rcursor, class);

    UNPROTECT(2);
    return Rcursor;
}


iproc_cursor *
Riproc_to_cursor (SEXP Rcursor)
{
    iproc_cursor *cursor;
    cursor = Riproc_sexp2ptr(Rcursor, FALSE, Riproc_cursor_type_tag, "cursor");
    return cursor;
}

SEXP
Riproc_cursor_new (SEXP Rmessages)
{
    iproc_messages *messages = Riproc_to_messages(Rmessages);
    iproc_cursor *cursor = iproc_cursor_new(messages);
    SEXP Rcursor;

    PROTECT(Rcursor = Riproc_from_cursor(cursor));
    iproc_cursor_unref(cursor);
    UNPROTECT(1);
    return Rcursor;
}

SEXP
Riproc_cursor_advance (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);

    if (iproc_cursor_next(cursor)) {
        return ScalarLogical(TRUE);
    } else {
        return ScalarLogical(FALSE);
    }
}

SEXP
Riproc_cursor_reset (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);
    iproc_cursor_reset(cursor);
    return NULL_USER_OBJECT;
}

SEXP
Riproc_cursor_time (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);

    if (!iproc_cursor_started(cursor))
        error("cursor has not started");
    if (iproc_cursor_finished(cursor))
        error("cursor has finished");

    int64_t time = iproc_cursor_time(cursor);
    return ScalarInteger(time);
}

SEXP
Riproc_cursor_nties (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);

    if (!iproc_cursor_started(cursor))
        error("cursor has not started");
    if (iproc_cursor_finished(cursor))
        error("cursor has finished");

    int64_t nties = iproc_cursor_nmsg(cursor);
    return ScalarInteger(nties);
}

SEXP
Riproc_cursor_from (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);

    if (!iproc_cursor_started(cursor))
        error("cursor has not started");
    if (iproc_cursor_finished(cursor))
        error("cursor has finished");

    int64_t i, n = iproc_cursor_nmsg(cursor);
    int *from;
    SEXP Rfrom;

    PROTECT(Rfrom = NEW_INTEGER(n));
    from = INTEGER_POINTER(Rfrom);

    for (i = 0; i < n; i++) {
        iproc_cursor_select_msg(cursor, i);
        int64_t msg_from = iproc_cursor_msg_from(cursor);
        from[i] = msg_from + 1;
    }

    UNPROTECT(1);
    return Rfrom;
}

SEXP
Riproc_cursor_to (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);

    if (!iproc_cursor_started(cursor))
        error("cursor has not started");
    if (iproc_cursor_finished(cursor))
        error("cursor has finished");

    int64_t i, n = iproc_cursor_nmsg(cursor);
    SEXP Rto;

    PROTECT(Rto = NEW_LIST(n));

    for (i = 0; i < n; i++) {
        iproc_cursor_select_msg(cursor, i);
        int64_t  msg_nto = iproc_cursor_msg_nto(cursor);
        int64_t *msg_to = iproc_cursor_msg_to(cursor);
        SEXP Rmsg_to;

        PROTECT(Rmsg_to = NEW_INTEGER(msg_nto));
        int64_t j;
        for (j = 0; j < msg_nto; j++) {
            INTEGER(Rmsg_to)[j] = msg_to[j] + 1;
        }
        SET_VECTOR_ELT(Rto, i, Rmsg_to);

        UNPROTECT(1);
    }

    UNPROTECT(1);
    return Rto;
}

SEXP
Riproc_cursor_started (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);
    if (iproc_cursor_started(cursor)) {
        return ScalarLogical(TRUE);
    } else {
        return ScalarLogical(FALSE);
    }
}

SEXP
Riproc_cursor_finished (SEXP Rcursor)
{
    iproc_cursor *cursor = Riproc_to_cursor(Rcursor);
    if (iproc_cursor_finished(cursor)) {
        return ScalarLogical(TRUE);
    } else {
        return ScalarLogical(FALSE);
    }
}
