#include "port.h"
#include <assert.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "r-utils.h"
#include "r-actors.h"
#include "r-cursor.h"
#include "r-design.h"
#include "r-model.h"

static SEXP Riproc_model_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_model_new", (DL_FUNC) & Riproc_model_new, 3},
	{"Riproc_model_design", (DL_FUNC) & Riproc_model_design, 1},
	{"Riproc_model_coefs", (DL_FUNC) & Riproc_model_coefs, 1},
	{"Riproc_model_has_loops", (DL_FUNC) & Riproc_model_has_loops, 1},
	{"Riproc_model_dim", (DL_FUNC) & Riproc_model_dim, 1},
	{"Riproc_model_nreceiver", (DL_FUNC) & Riproc_model_nreceiver, 1},
	{"Riproc_model_nsender", (DL_FUNC) & Riproc_model_nsender, 1},
	{"Riproc_model_log_probs", (DL_FUNC) & Riproc_model_log_probs, 3},
	{NULL, NULL, 0}
};

void Riproc_model_init(DllInfo * info)
{
	Riproc_model_type_tag = install("Riproc_model_type_tag");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void Riproc_model_free(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	iproc_model_unref(model);
}

iproc_model *Riproc_to_model(SEXP Rmodel)
{
	iproc_model *model =
	    Riproc_sexp2ptr(Rmodel, FALSE, Riproc_model_type_tag, "model");
	return model;
}

SEXP Riproc_from_model(iproc_model * model)
{
	SEXP Rmodel, class;

	iproc_model_ref(model);

	PROTECT(Rmodel =
		R_MakeExternalPtr(model, Riproc_model_type_tag, R_NilValue));
	R_RegisterCFinalizer(Rmodel, Riproc_model_free);

	/* set the class of the result */
	PROTECT(class = allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("model"));
	classgets(Rmodel, class);

	UNPROTECT(2);
	return Rmodel;
}

SEXP Riproc_model_new(SEXP Rdesign, SEXP Rcoefs, SEXP Rhas_loops)
{
	struct design *design = Riproc_to_design(Rdesign);
	struct vector coefs = Riproc_vector_view_sexp(Rcoefs);
	Rboolean has_loops = LOGICAL_VALUE(Rhas_loops);

	if (design_dim(design) != vector_dim(&coefs))
		error("design and coefs have different dimensions");
	if (has_loops == NA_LOGICAL)
		error("has.loops is be NA");

	iproc_model *model = iproc_model_new(design, &coefs, has_loops);
	SEXP Rmodel;

	PROTECT(Rmodel = Riproc_from_model(model));
	iproc_model_unref(model);

	UNPROTECT(1);
	return Rmodel;
}

SEXP Riproc_model_dim(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	int dim = (int)iproc_model_dim(model);
	return ScalarInteger(dim);
}

SEXP Riproc_model_nsender(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	int n = (int)iproc_model_nsender(model);
	return ScalarInteger(n);
}

SEXP Riproc_model_nreceiver(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	int n = (int)iproc_model_nreceiver(model);
	return ScalarInteger(n);
}

SEXP Riproc_model_design(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	struct design *design = iproc_model_design(model);
	return Riproc_from_design(design);
}

SEXP Riproc_model_coefs(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	struct vector *coefs = iproc_model_coefs(model);

	return Riproc_vector_new_copy(coefs);
}

SEXP Riproc_model_has_loops(SEXP Rmodel)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	if (iproc_model_has_loops(model)) {
		return ScalarLogical(TRUE);
	} else {
		return ScalarLogical(FALSE);
	}
}

SEXP Riproc_model_log_probs(SEXP Rmodel, SEXP Risend, SEXP Rcursor)
{
	iproc_model *model = Riproc_to_model(Rmodel);
	int i, n = GET_LENGTH(Risend);
	struct messages_iter *cursor = (Rcursor == NULL_USER_OBJECT
				      ? NULL : Riproc_to_cursor(Rcursor));
	int nsender = (int)iproc_model_nsender(model);
	int nreceiver = (int)iproc_model_nreceiver(model);
	struct history *history = messages_iter_history(cursor);

	SEXP Rprobst;

	PROTECT(Rprobst = allocMatrix(REALSXP, nreceiver, n));
	struct matrix probst = Riproc_matrix_view_sexp(Rprobst);
	struct vector dst;

	for (i = 0; i < n; i++) {
		int isend = INTEGER(Risend)[i] - 1;
		if (isend < 0 || isend >= nsender)
			error("invalid sender");

		iproc_model_ctx *ctx =
		    iproc_model_ctx_new(model, isend, history);
		vector_init_matrix_col(&dst, &probst, i);

		iproc_model_ctx_get_logprobs(ctx, &dst);
		iproc_model_ctx_unref(ctx);
	}

	UNPROTECT(1);
	return Rprobst;
}
