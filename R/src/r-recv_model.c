#include "port.h"

#include <stddef.h>
#include <errno.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "coreutil.h"
#include "history.h"


static int get_ids(size_t *dst, SEXP src)
{
	int i, n;
	int *xsrc;
	int err = 0;

	xsrc = INTEGER_POINTER(src);
	n = GET_LENGTH(src);

	for (i = 0; i < n; i++) {
		if (xsrc[i] <= 0) {
			err = -EDOM;
			goto out;
		}
		dst[i] = (size_t)(i - 1);
	}

out:
	return err;
}



SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver,
		       SEXP message_attr, SEXP nsend, SEXP nrecv,
		       SEXP loops, SEXP factors, SEXP types, SEXP traits,
		       SEXP specials)
{
	struct history h;
	double *xtime;
	double t;
	int *xsender, *xr;
	SEXP r;
	int i, n, j, m, xnsend, xnrecv, xloops;
	size_t from, *to;
	size_t nto, nto_max = 0;
	intptr_t attr;
	int nprotect = 0;

	n = GET_LENGTH(time);

	/* validate and extract 'time' */
	if (!IS_NUMERIC(time))
		error_return("'time' should be a numeric vector");

	xtime = NUMERIC_POINTER(time);

	for (i = 0; i < n; i++) {
		if (!R_FINITE(xtime[i])) {
			error_return("'time' values must be finite and non-missing");
		}
	}


	/* validate and extract 'nsend' */
	if (!IS_INTEGER(nsend) || GET_LENGTH(nsend) != 1)
		error_return("'nsend' should be a single integer");

	xnsend = INTEGER_VALUE(nsend);

	if (xnsend <= 0)
		error_return("'nsend' should be positive");


	/* validate and extract 'sender' */
	if (!IS_INTEGER(sender))
		error_return("'sender' should be an integer vector");
	if (GET_LENGTH(sender) != n)
		error_return("'time' and 'receiver' lengths differ");

	xsender = INTEGER_POINTER(sender);

	for (i = 0; i < n; i++) {
		if (!(1 <= xsender[i] && xsender[i] <= xnsend)) {
			error_return("'sender' value is out of range")
		}
	}


	/* validate and extract 'nrecv' */
	if (!IS_INTEGER(nrecv) || GET_LENGTH(nrecv) != 1)
		error_return("'nrecv' should be a single integer");

	xnrecv = INTEGER_VALUE(nrecv);

	if (xnrecv <= 0)
		error_return("'nrecv' should be positive");


	/* validate and extract 'receiver' */
	if (!IS_VECTOR(receiver))
		error_return("'receiver' should be a list");
	if (GET_LENGTH(receiver) != n)
		error_return("'time' and 'receiver' lengths differ");

	for (i = 0; i < n; i++) {
		r = VECTOR_ELT(receiver, i);
		if (!IS_INTEGER(r) || GET_LENGTH(r) == 0)
			error_return("each element of 'receiver' should be a non-empty integer vector");

		xr = INTEGER_POINTER(r);
		m = GET_LENGTH(r);

		nto_max = MAX(nto_max, (size_t)m);

		for (j = 0; j < m; j++) {
			if (!(1 <= xr[j] && xr[j] <= xnrecv)) {
				error_return("'receiver' value is out of range")
			}
		}
	}


	/* validate 'message.attr' */
	if (message_attr != NULL_USER_OBJECT)
		warning("'message.attr' has no effect");
	if (message_attr != NULL_USER_OBJECT && GET_LENGTH(message_attr) != n)
		error_return("'time' and 'message.attr' lengths differ");


	/* validate and extract 'loops' */
	if (!IS_LOGICAL(loops) || GET_LENGTH(loops) != 1)
		error_return("'loops' should be a single integer");

	xloops = LOGICAL_VALUE(loops);

	if (xloops == NA_LOGICAL)
		error_return("'loops' must be TRUE or FALSE")
	if (!xloops && xnrecv == 1)
		error_return("'nrecv' should be at least 2 (no loops)");


	to = (size_t *)R_alloc(nto_max, sizeof(size_t));

	/* create history object */
	history_init(&h, (size_t)xnsend, (size_t)xnrecv);
	for (i = 0; i < n; i++) {
		t = xtime[i];
		from = (size_t)(xsender[i] - 1);

		r = VECTOR_ELT(receiver, i);
		nto = GET_LENGTH(r);
		get_ids(to, r);

		attr = 0;

		history_set_time(&h, t);
		history_add(&h, from, to, nto, attr);
	}





	history_deinit(&h);
	error("Not implemented!");

	UNPROTECT(nprotect);
	return NULL_USER_OBJECT;
}

