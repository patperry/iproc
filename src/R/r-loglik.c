
#include "messages.h"
#include "r-messages.h"
#include "r-model.h"
#include "r-cursor.h"
#include "r-utils.h"
#include "r-loglik.h"


static SEXP Riproc_loglik_type_tag;

static R_CallMethodDef callMethods[] = {
    { "Riproc_loglik_new",    (DL_FUNC) &Riproc_loglik_new,    2 },
    { "Riproc_loglik_insert", (DL_FUNC) &Riproc_loglik_insert, 2 },
    { "Riproc_loglik_value",  (DL_FUNC) &Riproc_loglik_value,  1 },
    { "Riproc_loglik_grad",   (DL_FUNC) &Riproc_loglik_grad,   1 },
    { NULL,                   NULL,                            0 }
};


void
Riproc_loglik_init (DllInfo *info)
{
    Riproc_loglik_type_tag = install("Riproc_loglik_type_tag");
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void
Riproc_loglik_free (SEXP Rloglik)
{
    iproc_loglik *loglik = Riproc_to_loglik(Rloglik);
    iproc_loglik_unref(loglik);
}



iproc_loglik *
Riproc_to_loglik (SEXP Rloglik)
{
    iproc_loglik *loglik = Riproc_sexp2ptr(Rloglik, FALSE, Riproc_loglik_type_tag,
                                           "loglik");
    return loglik;
}


SEXP
Riproc_from_loglik (iproc_loglik *loglik)
{
    SEXP Rloglik, class;

    iproc_loglik_ref(loglik);

    PROTECT(Rloglik = R_MakeExternalPtr(loglik, Riproc_loglik_type_tag, R_NilValue));
    R_RegisterCFinalizer(Rloglik, Riproc_loglik_free);

    PROTECT(class = allocVector(STRSXP, 1));
    SET_STRING_ELT(class, 0, mkChar("loglik"));
    classgets(Rloglik, class);

    UNPROTECT(2);
    return Rloglik;
}


SEXP
Riproc_loglik_new (SEXP Rmodel,
                   SEXP Rmessages)
{
    iproc_model *model = Riproc_to_model(Rmodel);
    iproc_messages *messages = (Rmessages == NULL_USER_OBJECT
                                ? NULL
                                : Riproc_to_messages(Rmessages));
    iproc_loglik *loglik = iproc_loglik_new(model);
    SEXP Rloglik;
    PROTECT(Rloglik = Riproc_from_loglik(loglik));
    iproc_loglik_unref(loglik);

    if (messages) {
        iproc_message_iter *cursor = iproc_message_iter_new(messages);

        while (iproc_message_iter_next(cursor)) {
            iproc_history *history = iproc_message_iter_history(cursor);
            int64_t i, n = iproc_message_iter_ntie(cursor);

            for (i = 0; i < n; i++) {
                iproc_message_iter_select(cursor, i);
                int64_t  msg_from = iproc_message_iter_from(cursor);
                int64_t *msg_to = iproc_message_iter_to(cursor);
                int64_t  msg_nto = iproc_message_iter_nto(cursor);
                
                iproc_loglik_insertm(loglik, history, msg_from, msg_to, msg_nto);
            }
        }
        
        iproc_message_iter_unref(cursor);
    }

    UNPROTECT(1);
    return Rloglik;
}

SEXP
Riproc_loglik_insert (SEXP Rloglik,
                      SEXP Rcursor)
{
    iproc_loglik *loglik = Riproc_to_loglik(Rloglik);
    iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);
    iproc_history *history = iproc_message_iter_history(cursor);
    int64_t i, n = iproc_message_iter_ntie(cursor);

    for (i = 0; i < n; i++) {
        iproc_message_iter_select(cursor, i);
        int64_t  msg_from = iproc_message_iter_from(cursor);
        int64_t *msg_to = iproc_message_iter_to(cursor);
        int64_t  msg_nto = iproc_message_iter_nto(cursor);

        iproc_loglik_insertm(loglik, history, msg_from, msg_to, msg_nto);
    }

    return NULL_USER_OBJECT;
}

SEXP
Riproc_loglik_value (SEXP Rloglik)
{
    iproc_loglik *loglik = Riproc_to_loglik(Rloglik);
    double value = iproc_loglik_value(loglik);
    return ScalarReal(value);
}

SEXP
Riproc_loglik_grad  (SEXP Rloglik)
{
    iproc_loglik *loglik = Riproc_to_loglik(Rloglik);
    iproc_vector *grad = iproc_loglik_grad(loglik);
    return Riproc_vector_new_copy(grad);
}
