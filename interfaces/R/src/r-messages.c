#include "port.h"
#include <assert.h>
#include <stdint.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include "darray.h"
#include "r-utils.h"
#include "r-messages.h"

static SEXP Riproc_messages_type_tag;

static R_CallMethodDef callMethods[] = {
	{"Riproc_messages_new", (DL_FUNC) & Riproc_messages_new, 3},
	{"Riproc_messages_size", (DL_FUNC) & Riproc_messages_size, 1},
	{"Riproc_messages_time", (DL_FUNC) & Riproc_messages_time, 1},
	{"Riproc_messages_from", (DL_FUNC) & Riproc_messages_from, 1},
	{"Riproc_messages_to", (DL_FUNC) & Riproc_messages_to, 1},
	{"Riproc_messages_nto", (DL_FUNC) & Riproc_messages_nto, 1},
	{"Riproc_messages_max_from", (DL_FUNC) & Riproc_messages_max_from, 1},
	{"Riproc_messages_max_to", (DL_FUNC) & Riproc_messages_max_to, 1},
	{"Riproc_messages_max_nto", (DL_FUNC) & Riproc_messages_max_nto, 1},
	{NULL, NULL, 0}
};

void Riproc_messages_init(DllInfo * info)
{
	Riproc_messages_type_tag = install("Riproc_messages_type_tag");
	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}

static void Riproc_messages_free(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	iproc_messages_unref(msgs);
}

iproc_messages *Riproc_to_messages(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_sexp2ptr(Rmsgs,
					       FALSE,
					       Riproc_messages_type_tag,
					       "messages");
	return msgs;
}

SEXP Riproc_from_messages(iproc_messages * msgs)
{
	SEXP Rmsgs, class;

	iproc_messages_ref(msgs);

	PROTECT(Rmsgs =
		R_MakeExternalPtr(msgs, Riproc_messages_type_tag, R_NilValue));
	R_RegisterCFinalizer(Rmsgs, Riproc_messages_free);

	PROTECT(class = allocVector(STRSXP, 1));
	SET_STRING_ELT(class, 0, mkChar("messages"));
	classgets(Rmsgs, class);

	UNPROTECT(2);
	return Rmsgs;
}

static int64_t *copy_sexp_to_int64(struct darray *dst, SEXP Rsrc)
{
	int i, n = GET_LENGTH(Rsrc);
	int *src = INTEGER_POINTER(Rsrc);
	int s;

	darray_resize(dst, n);

	for (i = 0; i < n; i++) {
		s = src[i];
		if (s <= 0)
			error("'to' values must be positive");

		*(int64_t *)darray_at(dst, i) = s - 1;
	}

	if (n > 0) {
		return darray_front(dst);
	} else {
		return NULL;
	}
}

static void from_array_copy(int *dst, int64_t *src, int64_t n)
{
	int64_t i;
	int64_t s;

	for (i = 0; i < n; i++) {
		s = src[i];
		dst[i] = (int)s + 1;
	}
}

SEXP Riproc_messages_new(SEXP Rtime, SEXP Rfrom, SEXP Rto)
{
	int n = GET_LENGTH(Rtime);
	double *time = NUMERIC_POINTER(Rtime);
	int *from = INTEGER_POINTER(Rfrom);

	if (!(GET_LENGTH(Rfrom) == n && GET_LENGTH(Rto) == n))
		error("'time', 'from', and 'to' do not have same lengths");

	struct darray to_buf;
	double tcur = -INFINITY;
	double msg_time;
	int msg_from, msg_nto;
	int64_t *msg_to;
	int i;
	iproc_messages *msgs = iproc_messages_new();
	SEXP Rmsgs, Rmsg_to;

	darray_init(&to_buf, sizeof(int64_t));

	for (i = 0; i < n; i++) {
		msg_time = time[i];
		msg_from = from[i] - 1;
		Rmsg_to = VECTOR_ELT(Rto, i);
		msg_to = copy_sexp_to_int64(&to_buf, Rmsg_to);
		msg_nto = (int)darray_size(&to_buf);

		if (msg_time < tcur)
			error
			    ("'time' values must be sorted in increasing order");
		if (msg_from < 0)
			error("'from' values must be positive");

		iproc_messages_advance_to(msgs, msg_time);
		iproc_messages_insertm(msgs, msg_from, msg_to, msg_nto);

		tcur = msg_time;
	}

	PROTECT(Rmsgs = Riproc_from_messages(msgs));
	iproc_messages_unref(msgs);
	darray_deinit(&to_buf);

	UNPROTECT(1);
	return Rmsgs;
}

SEXP Riproc_messages_size(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int size = (int)iproc_messages_size(msgs);
	return ScalarInteger(size);
}

SEXP Riproc_messages_time(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int i, ntie, n = (int)iproc_messages_size(msgs);
	double t;
	SEXP Rtime;

	PROTECT(Rtime = NEW_NUMERIC(n));
	double *time = NUMERIC_POINTER(Rtime);

	iproc_message_iter *it = iproc_message_iter_new(msgs);

	while (iproc_message_iter_next(it)) {
		t = iproc_message_iter_time(it);
		ntie = (int)iproc_message_iter_ntie(it);
		for (i = 0; i < ntie; i++) {
			*time++ = t;
		}
	}

	iproc_message_iter_unref(it);

	UNPROTECT(1);
	return Rtime;
}

SEXP Riproc_messages_from(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int i, ntie, n = (int)iproc_messages_size(msgs);
	SEXP Rfrom;

	PROTECT(Rfrom = NEW_INTEGER(n));
	int *from = INTEGER_POINTER(Rfrom);

	iproc_message_iter *it = iproc_message_iter_new(msgs);

	while (iproc_message_iter_next(it)) {
		ntie = (int)iproc_message_iter_ntie(it);
		for (i = 0; i < ntie; i++) {
			iproc_message_iter_select(it, i);
			*from++ = (int)iproc_message_iter_from(it) + 1;
		}
	}

	iproc_message_iter_unref(it);

	UNPROTECT(1);
	return Rfrom;
}

SEXP Riproc_messages_to(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int i, msg, ntie, nto, n = (int)iproc_messages_size(msgs);
	int64_t *msg_to;
	SEXP Rto, Rmsg_to;
	iproc_message_iter *it;

	PROTECT(Rto = NEW_LIST(n));

	it = iproc_message_iter_new(msgs);
	msg = 0;

	while (iproc_message_iter_next(it)) {
		ntie = (int)iproc_message_iter_ntie(it);
		for (i = 0; i < ntie; i++) {
			iproc_message_iter_select(it, i);
			nto = (int)iproc_message_iter_nto(it);
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

SEXP Riproc_messages_nto(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int i, ntie, n = (int)iproc_messages_size(msgs);
	SEXP Rnto;

	PROTECT(Rnto = NEW_INTEGER(n));
	int *nto = INTEGER_POINTER(Rnto);

	iproc_message_iter *it = iproc_message_iter_new(msgs);

	while (iproc_message_iter_next(it)) {
		ntie = (int)iproc_message_iter_ntie(it);
		for (i = 0; i < ntie; i++) {
			iproc_message_iter_select(it, i);
			*nto++ = (int)iproc_message_iter_nto(it);
		}
	}

	iproc_message_iter_unref(it);

	UNPROTECT(1);
	return Rnto;
}

SEXP Riproc_messages_max_from(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int max_from = (int)iproc_messages_max_from(msgs);
	return ScalarInteger(max_from + 1);
}

SEXP Riproc_messages_max_to(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int max_to = (int)iproc_messages_max_to(msgs);
	return ScalarInteger(max_to + 1);
}

SEXP Riproc_messages_max_nto(SEXP Rmsgs)
{
	iproc_messages *msgs = Riproc_to_messages(Rmsgs);
	int max_nto = (int)iproc_messages_max_nto(msgs);
	return ScalarInteger(max_nto);
}
