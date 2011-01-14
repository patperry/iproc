#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <assert.h>
#include <stdint.h>
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

static void
Riproc_messages_free (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    iproc_messages_unref(msgs);
}


iproc_messages *
Riproc_to_messages (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_sexp2ptr(Rmsgs,
                                           FALSE,
                                           Riproc_messages_type_tag,
                                           "messages");
    return msgs;
}


SEXP
Riproc_from_messages (iproc_messages *msgs)
{
    SEXP Rmsgs, class;

    iproc_messages_ref(msgs);

    PROTECT(Rmsgs = R_MakeExternalPtr(msgs, Riproc_messages_type_tag, R_NilValue));
    R_RegisterCFinalizer(Rmsgs, Riproc_messages_free);

    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("messages"));
    classgets(Rmsgs, class);

    UNPROTECT(2);
    return Rmsgs;
}


static int
int_array_max (int *data, int n)
{
    int i;
    int max = 0;
    
    for (i = 0; i < n; i++) {
        int x = data[i];
        if (x > max)
            max = x;
    }

    return max;
}

static void
to_array_copy (int64_t *dst, int *src, int64_t n)
{
    int64_t i;
    int s;

    for (i = 0; i < n; i++) {
        s = src[i];
        if (s <= 0)
            error("'to' values must be positive");

        dst[i] = s - 1;
    }
}

static void
from_array_copy (int *dst, int64_t *src, int64_t n)
{
    int64_t i;
    int64_t s;

    for (i = 0; i < n; i++) {
        s = src[i];
        dst[i] = s + 1;
    }
}

SEXP
Riproc_messages_new (SEXP Rtime,
                     SEXP Rfrom,
                     SEXP Rto,
                     SEXP Rnto)
{
    int64_t n = GET_LENGTH(Rtime);
    int *time = INTEGER_POINTER(Rtime);
    int *from = INTEGER_POINTER(Rfrom);
    int *nto = INTEGER_POINTER(Rnto);
    int max_nto = int_array_max(nto, n);

    if (!(GET_LENGTH(Rfrom) == n && GET_LENGTH(Rnto) == n))
        error("'time', 'from', and 'nto' do not have same lengths");
    if (max_nto < 0)
        error("'nto' values must be nonnegative");

    int size_to = GET_LENGTH(Rto);
    int *to = INTEGER_POINTER(Rto);
    int *to_end = to + size_to;
    int64_t *msg_to = (int64_t *)R_alloc(max_nto, sizeof(int64_t));
    int64_t tcur = INT64_MIN;
    int64_t msg_from, msg_time, msg_nto;
    int64_t i;
    iproc_messages *msgs = iproc_messages_new(tcur);
    SEXP Rmsgs;

    for (i = 0; i < n; i++) {
        msg_time = time[i];
        msg_from = from[i] - 1;
        msg_nto = nto[i];

        if (msg_time < tcur)
            error("'time' values must be sorted in increasing order");
        if (msg_from < 0)
            error("'from' values must be positive");
        if (msg_nto < 0)
            error("'nto' values must be nonnegative");
        if (to + msg_nto > to_end)
            error("not enough elements in 'to' array");

        to_array_copy(msg_to, to, msg_nto);
        iproc_messages_advance_to(msgs, msg_time);
        iproc_messages_insertm(msgs, msg_from, msg_to, msg_nto);
        
        tcur = msg_time;
        to += msg_nto;
    }

    if (to != to_end)
        error("too many elements in 'to' array");
    
    PROTECT(Rmsgs = Riproc_from_messages(msgs));
    iproc_messages_unref(msgs);

    UNPROTECT(1);
    return Rmsgs;
}


SEXP
Riproc_messages_size (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t size = iproc_messages_size(msgs);
    return ScalarInteger(size);
}

SEXP
Riproc_messages_time (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t i, ntie, n = iproc_messages_size(msgs);
    int t;
    SEXP Rtime;

    PROTECT(Rtime = NEW_INTEGER(n));
    int *time = INTEGER_POINTER(Rtime);
    
    iproc_message_iter *it = iproc_message_iter_new(msgs);

    while (iproc_message_iter_next(it)) {
        t = iproc_message_iter_time(it);
        ntie = iproc_message_iter_ntie(it);
        for (i = 0; i < ntie; i++) {
            *time++ = t;
        }
    }

    iproc_message_iter_unref(it);

    UNPROTECT(1);
    return Rtime;
}

SEXP
Riproc_messages_from (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t i, ntie, n = iproc_messages_size(msgs);
    SEXP Rfrom;

    PROTECT(Rfrom = NEW_INTEGER(n));
    int *from = INTEGER_POINTER(Rfrom);
    
    iproc_message_iter *it = iproc_message_iter_new(msgs);

    while (iproc_message_iter_next(it)) {
        ntie = iproc_message_iter_ntie(it);
        for (i = 0; i < ntie; i++) {
            *from++ = iproc_message_iter_from(it) + 1;
        }
    }

    iproc_message_iter_unref(it);

    UNPROTECT(1);
    return Rfrom;
}

SEXP
Riproc_messages_to (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t i, msg, ntie, nto, n = iproc_messages_size(msgs);
    int64_t *msg_to;
    SEXP Rto, Rmsg_to;
    iproc_message_iter *it;

    PROTECT(Rto = NEW_LIST(n));

    it = iproc_message_iter_new(msgs);
    msg = 0;

    while (iproc_message_iter_next(it)) {
        ntie = iproc_message_iter_ntie(it);
        for (i = 0; i < ntie; i++) {
            nto = iproc_message_iter_nto(it);
            msg_to = iproc_message_iter_to(it);

            PROTECT(Rmsg_to = NEW_INTEGER(nto));
            from_array_copy(INTEGER_POINTER(Rmsg_to), msg_to, nto);
            SET_VECTOR_ELT(Rto, msg, Rmsg_to);
            UNPROTECT(1);
            msg++;
        }
    }
    assert(msg == n);

    iproc_message_iter_unref(it);

    UNPROTECT(1);
    return Rto;
}

SEXP
Riproc_messages_nto (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t i, ntie, n = iproc_messages_size(msgs);
    SEXP Rnto;

    PROTECT(Rnto = NEW_INTEGER(n));
    int *nto = INTEGER_POINTER(Rnto);
    
    iproc_message_iter *it = iproc_message_iter_new(msgs);

    while (iproc_message_iter_next(it)) {
        ntie = iproc_message_iter_ntie(it);
        for (i = 0; i < ntie; i++) {
            *nto++ = iproc_message_iter_nto(it);
        }
    }

    iproc_message_iter_unref(it);

    UNPROTECT(1);
    return Rnto;
}

SEXP
Riproc_messages_max_from (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t max_from = iproc_messages_max_from(msgs);
    return ScalarInteger(max_from + 1);
}
    
SEXP
Riproc_messages_max_to (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t max_to = iproc_messages_max_to(msgs);
    return ScalarInteger(max_to + 1);
}

SEXP
Riproc_messages_max_nto (SEXP Rmsgs)
{
    iproc_messages *msgs = Riproc_to_messages(Rmsgs);
    int64_t max_nto = iproc_messages_max_nto(msgs);
    return ScalarInteger(max_nto);
}
