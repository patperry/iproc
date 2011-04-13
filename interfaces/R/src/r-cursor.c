#include "port.h"
#include <assert.h>

#include "r-messages.h"
#include "r-utils.h"
#include "r-cursor.h"

static SEXP Riproc_cursor_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_cursor_new", (DL_FUNC) & Riproc_cursor_new, 2},
	{"Riproc_cursor_advance", (DL_FUNC) & Riproc_cursor_advance, 1},
	{"Riproc_cursor_started", (DL_FUNC) & Riproc_cursor_started, 1},
	{"Riproc_cursor_finished", (DL_FUNC) & Riproc_cursor_finished, 1},
	{"Riproc_cursor_reset", (DL_FUNC) & Riproc_cursor_reset, 1},
	{"Riproc_cursor_time", (DL_FUNC) & Riproc_cursor_time, 1},
	{"Riproc_cursor_nties", (DL_FUNC) & Riproc_cursor_nties, 1},
	{"Riproc_cursor_from", (DL_FUNC) & Riproc_cursor_from, 1},
	{"Riproc_cursor_to", (DL_FUNC) & Riproc_cursor_to, 1},
	{NULL, NULL, 0}
};

void Riproc_cursor_init(DllInfo * info)
{
	Riproc_cursor_type_tag = install("Riproc_cursor_type_tag");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void Riproc_cursor_free(SEXP Rcursor)
{
	iproc_message_iter *it = Riproc_to_cursor(Rcursor);
	iproc_message_iter_unref(it);
}

SEXP Riproc_from_cursor(iproc_message_iter * it)
{
	SEXP Rcursor, class;

	iproc_message_iter_ref(it);

	/* store the cursor pointer in an external pointer */
	PROTECT(Rcursor =
		R_MakeExternalPtr(it, Riproc_cursor_type_tag, R_NilValue));
	R_RegisterCFinalizer(Rcursor, Riproc_cursor_free);

	/* set the class of the result */
	PROTECT(class = allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("cursor"));
	classgets(Rcursor, class);

	UNPROTECT(2);
	return Rcursor;
}

iproc_message_iter *Riproc_to_cursor(SEXP Rcursor)
{
	iproc_message_iter *it;
	it = Riproc_sexp2ptr(Rcursor, FALSE, Riproc_cursor_type_tag, "cursor");
	return it;
}

SEXP Riproc_cursor_new(SEXP Rmessages)
{
	iproc_messages *messages = Riproc_to_messages(Rmessages);
	iproc_message_iter *cursor = iproc_message_iter_new(messages);
	SEXP Rcursor;

	PROTECT(Rcursor = Riproc_from_cursor(cursor));
	iproc_message_iter_unref(cursor);
	UNPROTECT(1);
	return Rcursor;
}

SEXP Riproc_cursor_advance(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);

	if (iproc_message_iter_next(cursor)) {
		return ScalarLogical(TRUE);
	} else {
		return ScalarLogical(FALSE);
	}
}

SEXP Riproc_cursor_reset(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);
	iproc_message_iter_reset(cursor);
	return NULL_USER_OBJECT;
}

SEXP Riproc_cursor_time(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);

	if (!iproc_message_iter_started(cursor))
		error("cursor has not started");
	if (iproc_message_iter_finished(cursor))
		error("cursor has finished");

	double time = iproc_message_iter_time(cursor);
	return ScalarReal(time);
}

SEXP Riproc_cursor_nties(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);

	if (!iproc_message_iter_started(cursor))
		error("cursor has not started");
	if (iproc_message_iter_finished(cursor))
		error("cursor has finished");

	int nties = (int)iproc_message_iter_ntie(cursor);
	return ScalarInteger(nties);
}

SEXP Riproc_cursor_from(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);

	if (!iproc_message_iter_started(cursor))
		error("cursor has not started");
	if (iproc_message_iter_finished(cursor))
		error("cursor has finished");

	int i, n = (int)iproc_message_iter_ntie(cursor);
	int *from;
	SEXP Rfrom;

	PROTECT(Rfrom = NEW_INTEGER(n));
	from = INTEGER_POINTER(Rfrom);

	for (i = 0; i < n; i++) {
		iproc_message_iter_select(cursor, i);
		int msg_from = (int)iproc_message_iter_from(cursor);
		from[i] = msg_from + 1;
	}

	UNPROTECT(1);
	return Rfrom;
}

SEXP Riproc_cursor_to(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);

	if (!iproc_message_iter_started(cursor))
		error("cursor has not started");
	if (iproc_message_iter_finished(cursor))
		error("cursor has finished");

	int i, n = (int)iproc_message_iter_ntie(cursor);
	SEXP Rto;

	PROTECT(Rto = NEW_LIST(n));

	for (i = 0; i < n; i++) {
		iproc_message_iter_select(cursor, i);
		int msg_nto = (int)iproc_message_iter_nto(cursor);
		int64_t *msg_to = iproc_message_iter_to(cursor);
		SEXP Rmsg_to;

		PROTECT(Rmsg_to = NEW_INTEGER(msg_nto));
		int j;
		for (j = 0; j < msg_nto; j++) {
			INTEGER(Rmsg_to)[j] = (int)msg_to[j] + 1;
		}
		SET_VECTOR_ELT(Rto, i, Rmsg_to);

		UNPROTECT(1);
	}

	UNPROTECT(1);
	return Rto;
}

SEXP Riproc_cursor_started(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);
	if (iproc_message_iter_started(cursor)) {
		return ScalarLogical(TRUE);
	} else {
		return ScalarLogical(FALSE);
	}
}

SEXP Riproc_cursor_finished(SEXP Rcursor)
{
	iproc_message_iter *cursor = Riproc_to_cursor(Rcursor);
	if (iproc_message_iter_finished(cursor)) {
		return ScalarLogical(TRUE);
	} else {
		return ScalarLogical(FALSE);
	}
}
