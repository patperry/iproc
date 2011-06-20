#include "port.h"
#include <assert.h>

#include "messages.h"
#include "r-messages.h"
#include "r-model.h"
#include "r-utils.h"
#include "r-loglik.h"

static SEXP Riproc_loglik_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_loglik_new", (DL_FUNC) & Riproc_loglik_new, 2},
	{"Riproc_loglik_value", (DL_FUNC) & Riproc_loglik_value, 1},
	{"Riproc_loglik_grad", (DL_FUNC) & Riproc_loglik_grad, 1},
	{NULL, NULL, 0}
};

void Riproc_loglik_init(DllInfo * info)
{
	Riproc_loglik_type_tag = install("Riproc_loglik_type_tag");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void Riproc_loglik_free(SEXP Rloglik)
{
	struct recv_loglik *loglik = Riproc_to_loglik(Rloglik);
	recv_loglik_free(loglik);
}

struct recv_loglik *Riproc_to_loglik(SEXP Rloglik)
{
	struct recv_loglik *loglik =
	    Riproc_sexp2ptr(Rloglik, FALSE, Riproc_loglik_type_tag,
			    "loglik");
	return loglik;
}

SEXP Riproc_from_loglik(struct recv_loglik * loglik)
{
	SEXP Rloglik, class;

	// iproc_loglik_ref(loglik);

	PROTECT(Rloglik =
		R_MakeExternalPtr(loglik, Riproc_loglik_type_tag, R_NilValue));
	R_RegisterCFinalizer(Rloglik, Riproc_loglik_free);

	PROTECT(class = allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("loglik"));
	classgets(Rloglik, class);

	UNPROTECT(2);
	return Rloglik;
}

SEXP Riproc_loglik_new(SEXP Rmodel, SEXP Rmessages)
{
	struct model *model = Riproc_to_model(Rmodel);
	struct messages *messages = (Rmessages == NULL_USER_OBJECT
				     ? NULL : Riproc_to_messages(Rmessages));
	struct recv_loglik *loglik = recv_loglik_alloc(model, messages);
	SEXP Rloglik;
	PROTECT(Rloglik = Riproc_from_loglik(loglik));
	recv_loglik_free(loglik);
	UNPROTECT(1);
	return Rloglik;
}

SEXP Riproc_loglik_value(SEXP Rloglik)
{
	struct recv_loglik *loglik = Riproc_to_loglik(Rloglik);
	double value = recv_loglik_avg_dev(loglik);
	return ScalarReal(value);
}

SEXP Riproc_loglik_grad(SEXP Rloglik)
{
	//struct recv_loglik *loglik = Riproc_to_loglik(Rloglik);
	//struct vector *grad = recv_loglik_grad(loglik);
	//return Riproc_vector_new_copy(grad);
	return NULL;
}
