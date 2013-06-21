#include "port.h"

#include <stddef.h>
#include <errno.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "coreutil.h"

#include "history.h"
#include "design.h"


#define DOMAIN_ERROR(msg) \
	do { \
		error(msg); \
		return -EDOM; \
	} while (0)



struct args {
	/* dimensions */
	size_t nsend;
	size_t nrecv;

	/* message data */
	double *time;
	int *sender;
	SEXP receiver;
	size_t nmsg;
	size_t nto_max;

	/* properties */
	int loops;
};



static int extract_args(SEXP time, SEXP sender, SEXP receiver,
			SEXP message_attr, SEXP nsend, SEXP nrecv,
			SEXP loops, SEXP factors, SEXP types, SEXP traits,
			SEXP specials, struct args *args);
static void setup_history(struct history *h, const struct args *args);
static void setup_recv_design(struct design *r, struct history *h,
			      const struct args *args);
static int get_ids(size_t *dst, SEXP src);



SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver,
		       SEXP message_attr, SEXP nsend, SEXP nrecv,
		       SEXP loops, SEXP factors, SEXP types, SEXP traits,
		       SEXP specials)
{
	struct history history;
	struct design recv;
	struct args args;
	int err = 0;

	err = extract_args(time, sender, receiver, message_attr, nsend, nrecv,
			   loops, factors, types, traits, specials, &args);

	if (err < 0)
		goto args_fail;

	setup_history(&history, &args);
	setup_recv_design(&recv, &history, &args);



	design_deinit(&recv);
	history_deinit(&history);
	error("Not implemented!");

args_fail:
	return NULL_USER_OBJECT;
}




static int extract_args(SEXP time, SEXP sender, SEXP receiver,
			SEXP message_attr, SEXP nsend, SEXP nrecv,
			SEXP loops, SEXP factors, SEXP types, SEXP traits,
			SEXP specials, struct args *args)
{
	double *xtime;
	int *xsender, *xr;
	SEXP r;
	int i, n, j, m, xnsend, xnrecv, xloops;
	size_t nto_max = 0;

	n = GET_LENGTH(time);

	/* validate and extract 'time' */
	if (!IS_NUMERIC(time))
		DOMAIN_ERROR("'time' should be a numeric vector");

	xtime = NUMERIC_POINTER(time);

	for (i = 0; i < n; i++) {
		if (!R_FINITE(xtime[i]))
			DOMAIN_ERROR("'time' values must be finite and non-missing");
	}


	/* validate and extract 'nsend' */
	if (!IS_INTEGER(nsend) || GET_LENGTH(nsend) != 1)
		DOMAIN_ERROR("'nsend' should be a single integer");

	xnsend = INTEGER_VALUE(nsend);

	if (xnsend <= 0)
		DOMAIN_ERROR("'nsend' should be positive");


	/* validate and extract 'sender' */
	if (!IS_INTEGER(sender))
		DOMAIN_ERROR("'sender' should be an integer vector");
	if (GET_LENGTH(sender) != n)
		DOMAIN_ERROR("'time' and 'receiver' lengths differ");

	xsender = INTEGER_POINTER(sender);

	for (i = 0; i < n; i++) {
		if (!(1 <= xsender[i] && xsender[i] <= xnsend)) {
			DOMAIN_ERROR("'sender' value is out of range");
		}
	}


	/* validate and extract 'nrecv' */
	if (!IS_INTEGER(nrecv) || GET_LENGTH(nrecv) != 1)
		DOMAIN_ERROR("'nrecv' should be a single integer");

	xnrecv = INTEGER_VALUE(nrecv);

	if (xnrecv <= 0)
		DOMAIN_ERROR("'nrecv' should be positive");


	/* validate and extract 'receiver' */
	if (!IS_VECTOR(receiver))
		DOMAIN_ERROR("'receiver' should be a list");
	if (GET_LENGTH(receiver) != n)
		DOMAIN_ERROR("'time' and 'receiver' lengths differ");


	for (i = 0; i < n; i++) {
		r = VECTOR_ELT(receiver, i);
		if (!IS_INTEGER(r) || GET_LENGTH(r) == 0)
			DOMAIN_ERROR("each element of 'receiver' should be a non-empty integer vector");

		xr = INTEGER_POINTER(r);
		m = GET_LENGTH(r);

		nto_max = MAX(nto_max, (size_t)m);

		for (j = 0; j < m; j++) {
			if (!(1 <= xr[j] && xr[j] <= xnrecv)) {
				DOMAIN_ERROR("'receiver' value is out of range");
			}
		}
	}


	/* validate 'message.attr' */
	if (message_attr != NULL_USER_OBJECT)
		warning("'message.attr' has no effect");
	if (message_attr != NULL_USER_OBJECT && GET_LENGTH(message_attr) != n)
		DOMAIN_ERROR("'time' and 'message.attr' lengths differ");


	/* validate and extract 'loops' */
	if (!IS_LOGICAL(loops) || GET_LENGTH(loops) != 1)
		DOMAIN_ERROR("'loops' should be a single integer");

	xloops = LOGICAL_VALUE(loops);

	if (xloops == NA_LOGICAL)
		DOMAIN_ERROR("'loops' must be TRUE or FALSE");
	if (!xloops && xnrecv == 1)
		DOMAIN_ERROR("'nrecv' should be at least 2 (no loops)");

	args->nsend = (size_t)xnsend;
	args->nrecv = (size_t)xnrecv;
	args->time = xtime;
	args->sender = xsender;
	args->receiver = receiver;
	args->nmsg = (size_t)n;
	args->nto_max = nto_max;
	args->loops = xloops;

	return 0;
}


static void setup_history(struct history *h, const struct args *args)
{
	size_t nsend;
	size_t nrecv;
	double t, *time;
	int *sender;
	SEXP r, receiver;
	size_t from, *to;
	size_t nto, nto_max;
	intptr_t attr;
	size_t i, n;

	nsend = args->nsend;
	nrecv = args->nrecv;

	n = args->nmsg;
	time = args->time;
	sender = args->sender;
	receiver = args->receiver;
	nto_max = args->nto_max;
	to = (size_t *)R_alloc(nto_max, sizeof(size_t));


	history_init(h, nsend, nrecv);
	for (i = 0; i < n; i++) {
		t = time[i];
		from = (size_t)(sender[i] - 1);
		r = VECTOR_ELT(receiver, i);
		nto = GET_LENGTH(r);
		get_ids(to, r);

		attr = 0;

		history_set_time(h, t);
		history_add(h, from, to, nto, attr);
	}
}


static void setup_recv_design(struct design *r, struct history *h,
			      const struct args *args)
{
	size_t nrecv = args->nrecv;
	design_init(r, h, nrecv);
}


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
		dst[i] = (size_t)(xsrc[i] - 1);
	}

out:
	return err;
}

