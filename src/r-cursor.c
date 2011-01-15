
#include "r-messages.h"
#include "r-utils.h"
#include "r-cursor.h"



static SEXP Riproc_cursor_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_cursor_new",     (DL_FUNC) &Riproc_cursor_new,     2 },
    { "Riproc_cursor_advance", (DL_FUNC) &Riproc_cursor_advance, 1 },
    { "Riproc_cursor_reset",   (DL_FUNC) &Riproc_cursor_reset,   1 },
    { NULL,                    NULL,                             0 }
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
