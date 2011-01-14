#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-messages.h"

static SEXP Riproc_messages_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_messages_new",      (DL_FUNC) &Riproc_messages_new,      4 },
    { "Riproc_messages_size",     (DL_FUNC) &Riproc_messages_size,     1 },
    { "Riproc_messages_time",     (DL_FUNC) &Riproc_messages_time,     1 },
    { "Riproc_messages_from",     (DL_FUNC) &Riproc_messages_from,     1 },
    { "Riproc_messages_to",       (DL_FUNC) &Riproc_messages_to,       1 },
    { "Riproc_messages_nto",      (DL_FUNC) &Riproc_messages_nto,      1 },
    { "Riproc_messages_max_from", (DL_FUNC) &Riproc_messages_max_from, 1 },
    { "Riproc_messages_max_to",   (DL_FUNC) &Riproc_messages_max_to,   1 },
    { "Riproc_messages_max_nto",  (DL_FUNC) &Riproc_messages_max_nto,  1 },
    { NULL,                       NULL,                                0 }
};


void
Riproc_messages_init (DllInfo *info)
{
    Riproc_messages_type_tag = install("Riproc_messages_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

iproc_messages *
Riproc_to_messages (SEXP Rmsgs)
{

}

SEXP
Riproc_from_messages (iproc_messages *msgs)
{
}


SEXP
Riproc_messages_new (SEXP Rtime,
                     SEXP Rfrom,
                     SEXP Rto,
                     SEXP Rnto)
{
}

SEXP
Riproc_messages_size (SEXP Rmsgs)
{

}

SEXP
Riproc_messages_time (SEXP Rmsgs)
{

}

SEXP
Riproc_messages_from (SEXP Rmsgs)
{
}

SEXP
Riproc_messages_to (SEXP Rmsgs)
{
}

SEXP
Riproc_messages_nto (SEXP Rmsgs)
{
}

SEXP
Riproc_messages_max_from (SEXP Rmsgs)
{
}
    
SEXP
Riproc_messages_max_to (SEXP Rmsgs)
{
}

SEXP
Riproc_messages_max_nto (SEXP Rmsgs)
{
}


