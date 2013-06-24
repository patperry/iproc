#include "port.h"

#include <errno.h>
#include <stddef.h>
#include <string.h>

#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>

#include "blas.h"
#include "coreutil.h"
#include "matrixutil.h"

#include "history.h"
#include "design.h"
#include "design2.h"
#include "recv_fit.h"


#define CHECK(x) \
	do { \
		if ((err = x)) \
			return err; \
	} while (0)


#define DOMAIN_ERROR(msg) \
	do { \
		error(msg); \
		return -EDOM; \
	} while (0)


enum variable_type {
	VARIABLE_TYPE_TRAIT,
	VARIABLE_TYPE_SPECIAL,
	VARIABLE_TYPE_RESPONSE
};



struct args {
	/* dimensions */
	size_t nsend;
	size_t nrecv;

	/* properties */
	int loops;

	/* message data */
	double *time;
	int *sender;
	SEXP receiver; /* list of integer vectors */
	size_t nmsg;
	size_t nto_max;

	/* terms */
	size_t nvar;
	size_t nterm;
	enum variable_type *types;
	SEXP names; /* character vector */
	SEXP term_names; /* character vector */
	int *factors;
	size_t *trait_ind;
	size_t ntrait;
	size_t *special_ind;
	size_t nspecial;

	/* traits */
	SEXP traits;
	SEXP assign;
};



static int extract_args(SEXP time, SEXP sender, SEXP receiver,
			SEXP message_attr, SEXP nsend, SEXP nrecv,
			SEXP loops, SEXP factors, SEXP types, SEXP traits,
			SEXP specials, struct args *args);
static int extract_properties(SEXP nsend, SEXP nrecv, SEXP loops,
			      struct args *args);
static int extract_messages(SEXP time, SEXP sender, SEXP receiver,
			    SEXP message_attr, size_t nsend, size_t nrecv,
			    struct args *args);
static int extract_terms(SEXP factors, SEXP types, struct args *args);
static int extract_traits(SEXP traits, size_t nrecv, size_t ntrait, struct args *args);



static void setup_history(struct history *h, const struct args *args);
static void setup_recv_design(struct design *r, struct history *h,
			      const struct args *args);
static void setup_dyad_design(struct design2 *d, struct design *r,
			      struct history *h, const struct args *args);
static int do_fit(struct recv_fit *fit, const struct recv_params *params0, const double *duals0);


static int get_ids(size_t *dst, SEXP src);



SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver,
		       SEXP message_attr, SEXP nsend, SEXP nrecv,
		       SEXP loops, SEXP factors, SEXP types, SEXP traits,
		       SEXP specials)
{
	struct history history;
	struct design recv;
	struct design2 dyad;
	int exclude_loops;
	const struct message *msgs;
	size_t nmsg;
	size_t ncextra;
	struct recv_fit fit;
	struct args args;
	int err = 0;

	err = extract_args(time, sender, receiver, message_attr, nsend, nrecv,
			   loops, factors, types, traits, specials, &args);

	if (err < 0)
		goto args_fail;

	setup_history(&history, &args);
	setup_recv_design(&recv, &history, &args);
	setup_dyad_design(&dyad, &recv, &history, &args);

	history_get_messages(&history, &msgs, &nmsg);
	exclude_loops = args.loops ? 0 : 1;

	recv_fit_init(&fit, &recv, &dyad, exclude_loops, msgs, nmsg, NULL, NULL);

	ncextra = recv_fit_extra_constr_count(&fit);
	if (ncextra)
		warning("Adding %zd %s to make parameters identifiable\n",
			ncextra, ncextra == 1 ? "constraint" : "constraints");

	err = do_fit(&fit, NULL, NULL);


	recv_fit_deinit(&fit);
	design2_deinit(&dyad);
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
	int err;

	CHECK(extract_properties(nsend, nrecv, loops, args));
	CHECK(extract_messages(time, sender, receiver, message_attr,
			       args->nsend, args->nrecv, args));
	CHECK(extract_terms(factors, types, args));
	CHECK(extract_traits(traits, args->nrecv, args->ntrait, args));

	return 0;
}


static int extract_properties(SEXP nsend, SEXP nrecv, SEXP loops,
			      struct args *args)
{
	int xnsend, xnrecv, xloops;

	/* validate and extract 'nsend' */
	if (!IS_INTEGER(nsend) || GET_LENGTH(nsend) != 1)
		DOMAIN_ERROR("'nsend' should be a single integer");

	xnsend = INTEGER_VALUE(nsend);

	if (xnsend <= 0)
		DOMAIN_ERROR("'nsend' should be positive");


	/* validate and extract 'nrecv' */
	if (!IS_INTEGER(nrecv) || GET_LENGTH(nrecv) != 1)
		DOMAIN_ERROR("'nrecv' should be a single integer");

	xnrecv = INTEGER_VALUE(nrecv);

	if (xnrecv <= 0)
		DOMAIN_ERROR("'nrecv' should be positive");



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
	args->loops = xloops;

	return 0;
}


static int extract_messages(SEXP time, SEXP sender, SEXP receiver,
			    SEXP message_attr, size_t nsend, size_t nrecv,
			    struct args *args)
{
	double *xtime;
	int *xsender, *xr;
	SEXP r;
	int i, n, j, m;
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


	/* validate and extract 'sender' */
	if (!IS_INTEGER(sender))
		DOMAIN_ERROR("'sender' should be an integer vector");
	if (GET_LENGTH(sender) != n)
		DOMAIN_ERROR("'time' and 'receiver' lengths differ");

	xsender = INTEGER_POINTER(sender);

	for (i = 0; i < n; i++) {
		if ((size_t)xsender[i] > nsend) {
			DOMAIN_ERROR("'sender' value is out of range");
		}
	}



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
			if ((size_t)xr[j] > nrecv) {
				DOMAIN_ERROR("'receiver' value is out of range");
			}
		}
	}


	/* validate 'message.attr' */
	if (message_attr != NULL_USER_OBJECT)
		warning("'message.attr' has no effect");
	if (message_attr != NULL_USER_OBJECT && GET_LENGTH(message_attr) != n)
		DOMAIN_ERROR("'time' and 'message.attr' lengths differ");


	args->time = xtime;
	args->sender = xsender;
	args->receiver = receiver;
	args->nmsg = (size_t)n;
	args->nto_max = nto_max;

	return 0;
}


static int extract_terms(SEXP factors, SEXP types, struct args *args)
{
	SEXP dim, rl, cl;
	const char *rn, *cn, *xtype;
	int i, j, m, n;
	int *xfactors;
	enum variable_type *xtypes;
	size_t *trait_ind, *special_ind;
	size_t ntrait, nspecial;

	/* validate and extract factor matrix */
	if (!IS_INTEGER(factors))
		DOMAIN_ERROR("'factors' should be integer");
	if (!isMatrix(factors))
		DOMAIN_ERROR("'factors' should be a matrix");

	dim = GET_DIM(factors);
	m = INTEGER(dim)[0];
	n = INTEGER(dim)[1];
	xfactors = INTEGER_POINTER(factors);

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			if (xfactors[i * n + j] != 0 && xfactors[i * n + j] != 1)
				DOMAIN_ERROR("invalid entry in 'factors' matrix");
		}
	}


	/* validate and extract variable and term names */
	GetMatrixDimnames(factors, &rl, &cl, &rn, &cn);

	if (!IS_CHARACTER(rl))
		DOMAIN_ERROR("'factors' should have character rownames");
	if (GET_LENGTH(rl) != m)
		DOMAIN_ERROR("'factors' should have a name for each row");
	if (n > 0) {
		if (!IS_CHARACTER(cl))
			DOMAIN_ERROR("'factors' should have character colnames");
		if (GET_LENGTH(cl) != n)
			DOMAIN_ERROR("'factors' should have a name for each column");
	}


	/* validate and extract variable types */
	if (!IS_CHARACTER(types))
		DOMAIN_ERROR("'types' should be a character vector");
	if (LENGTH(types) != m)
		DOMAIN_ERROR("'types' vector has the wrong length");

	xtypes = (void *)R_alloc(m, sizeof(*xtypes));

	trait_ind = (void *)R_alloc(m, sizeof(*trait_ind));
	special_ind = (void *)R_alloc(m, sizeof(*special_ind));
	ntrait = 0;
	nspecial = 0;

	for (i = 0; i < m; i++) {
		xtype = CHAR(STRING_ELT(types, i));

		if (strcmp(xtype, "response") == 0) {
			xtypes[i] = VARIABLE_TYPE_RESPONSE;
		} else if (strcmp(xtype, "special") == 0) {
			xtypes[i] = VARIABLE_TYPE_SPECIAL;
			special_ind[nspecial++] = (size_t)i;
		} else if (strcmp(xtype, "trait") == 0) {
			xtypes[i] = VARIABLE_TYPE_TRAIT;
			trait_ind[ntrait++] = (size_t)i;
		} else {
			DOMAIN_ERROR("invalid 'types' value");
		}
	}

	args->nvar = m;
	args->nterm = n;
	args->types = xtypes;
	args->names = rl;
	args->term_names = cl;
	args->factors = xfactors;
	args->ntrait = ntrait;
	args->trait_ind = trait_ind;
	args->nspecial = nspecial;
	args->special_ind = special_ind;

	return 0;
}


static int extract_traits(SEXP traits, size_t nrecv, size_t ntrait, struct args *args)
{
	SEXP dim, assign;
	int i, m, n;
	int *xassign;

	/* validate and extract factor matrix */
	if (!IS_NUMERIC(traits))
		DOMAIN_ERROR("'traits' should be numeric");
	if (!isMatrix(traits))
		DOMAIN_ERROR("'traits' should be a matrix");

	dim = GET_DIM(traits);
	m = INTEGER(dim)[0];
	n = INTEGER(dim)[1];

	if ((size_t)m != nrecv)
		DOMAIN_ERROR("'traits' should have one row for each receiver");

	assign = GET_ATTR(traits, install("assign"));

	if (assign == NULL_USER_OBJECT)
		DOMAIN_ERROR("'traits' should have an 'assign' attribute");
	if (!IS_INTEGER(assign))
		DOMAIN_ERROR("'assign' should be integer");
	if (GET_LENGTH(assign) != n)
		DOMAIN_ERROR("'assign' should have length"
			     " equal to the number of columns in 'traits'");

	xassign = INTEGER_POINTER(assign);
	for (i = 0; i < n; i++) {
		if (!(0 <= xassign[i] && xassign[i] <= (int)ntrait))
			DOMAIN_ERROR("'assign' values should be between 0"
				     " and the number of trait variables");
	}

	args->traits = traits;
	args->assign = assign;

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
	SEXP traits_dims;
	double *traits;
	size_t nrecv;
	size_t i, ntrait;
	size_t j, dim_max;
	size_t dim, rank;
	size_t *dims;
	int *assign;
	const char *name = NULL;
	double *buf, *x;

	nrecv = args->nrecv;
	ntrait = args->ntrait;
	traits_dims = GET_DIM(args->traits);
	dim_max = (size_t)(INTEGER(traits_dims)[1]);

	traits = NUMERIC_POINTER(args->traits);
	assign = INTEGER_POINTER(args->assign);

	design_init(r, h, nrecv);

	buf = (void *)R_alloc(dim_max * nrecv, sizeof(*buf));
	x = (void *)R_alloc(nrecv * dim_max, sizeof(*x));

	for (i = 0; i < ntrait; i++) {
		dim = 0;
		for (j = 0; j < dim_max; j++) {
			if ((size_t)assign[j] == i + 1) {
				blas_dcopy(nrecv, traits + j * nrecv, 1,
					   buf + dim * nrecv, 1);
				dim++;
			}
		}

		if (dim == 0) {
			/* TODO: check for in extract_traits */
			error("0-dimensional trait %zd", i);
		}

		matrix_dtrans(nrecv, dim, buf, nrecv, x, dim);

		if (dim > 1) {
			rank = 1;
			dims = &dim;
		} else {
			rank = 0;
			dims = NULL;
		}

		design_add_trait(r, name, x, dims, rank);
	}
}


static void setup_dyad_design(struct design2 *d, struct design *r,
			      struct history *h, const struct args *args)
{
	size_t nsend, nrecv;

	nsend = args->nsend;
	nrecv = args->nrecv;

	design2_init(d, h, nsend, nrecv);
}



static int do_fit(struct recv_fit *fit, const struct recv_params *params0, const double *duals0)
{
	size_t maxit = 300;
	size_t report = 1;
	int trace = 1;
	enum recv_fit_task task;
	size_t it = 0;
	int err = 0;

	task = recv_fit_start(fit, params0, duals0);
	for (it = 0; it < maxit && task == RECV_FIT_STEP; it++) {
		task = recv_fit_advance(fit);

		if (trace && it % report == 0) {
			// const struct recv_loglik *ll = recv_fit_loglik(&fit);
			// const struct recv_loglik_info *info = recv_loglik_info(ll);
			// size_t n = info->nrecv;
			double dev = recv_fit_dev(fit);
			double nscore = recv_fit_score_norm(fit);
			double step = recv_fit_step_size(fit);

			Rprintf("iter %zu deviance = %.2f; |score| = %.16f; step = %.16f\n",
				it, dev, nscore, step);
		}
	}

	if (task != RECV_FIT_CONV) {
		error("%s", recv_fit_errmsg(task));
		err = -1;
	}

	return err;
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

