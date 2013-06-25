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
	VARIABLE_TYPE_RECV_TRAIT,
	VARIABLE_TYPE_RECV_SPECIAL,
	VARIABLE_TYPE_DYAD_SPECIAL,
	VARIABLE_TYPE_RESPONSE
};


struct variable {
	enum variable_type type;
	const char *name;
	const char **names;
	size_t size;
	union { struct var *recv; struct var2 *dyad; } var;
};

/*
enum term_type {
	TERM_TYPE_RECV_TRAIT,
	TERM_TYPE_RECV_SPECIAL,
	TERM_TYPE_RECV_PROD,
	TERM_TYPE_DYAD_SPECIAL,
};
*/


struct properties {
	/* dimensions */
	size_t nsend;
	size_t nrecv;

	/* properties */
	int exclude_loops;
};


struct terms {
	size_t nvar;
	size_t nterm;
};


struct args {


	/* terms */
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



static int get_properties(SEXP nsend, SEXP nrecv, SEXP loops,
			  struct properties *p);

static int get_history(SEXP time, SEXP sender, SEXP receiver,
		       struct history *h);

static int get_terms(SEXP factors, SEXP variables, SEXP order,
		     struct design *s, struct design *r, struct design2 *d,
		     struct terms *tm);

//static int extract_terms(SEXP factors, SEXP types, struct args *args);
//static int extract_traits(SEXP traits, size_t nrecv, size_t ntrait, struct args *args);



//static void setup_recv_design(struct design *r, struct history *h,
//			      const struct args *args);
//static void setup_dyad_design(struct design2 *d, struct design *r,
//			      struct history *h, const struct args *args);
//static int do_fit(struct recv_fit *fit, const struct recv_params *params0, const double *duals0);


//static int get_ids(size_t *dst, SEXP src);



SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver,
		       SEXP factors, SEXP variables, SEXP order,
		       SEXP nsend, SEXP nrecv, SEXP loops)
{
	struct properties p;
	struct history h;
	struct design s, r;
	struct design2 d;
	struct terms tm;

	const struct message *msgs;
	size_t nmsg;
	//size_t ncextra;
	//struct recv_fit fit;
	int err = 0;

	err = get_properties(nsend, nrecv, loops, &p);
	if (err < 0)
		goto properties_fail;

	history_init(&h, p.nsend, p.nrecv);
	err = get_history(time, sender, receiver, &h);
	if (err < 0)
		goto history_fail;

	design_init(&s, &h, p.nsend);
	design_init(&r, &h, p.nrecv);
	design2_init(&d, &h, p.nsend, p.nrecv);

	err = get_terms(factors, variables, order, &s, &r, &d, &tm);
	if (err < 0)
		goto terms_fail;

	//setup_recv_design(&recv, &history, &args);
	//setup_dyad_design(&dyad, &recv, &history, &args);

	history_get_messages(&h, &msgs, &nmsg);

	//recv_fit_init(&fit, &recv, &dyad, p.exclude_loops, msgs, nmsg, NULL, NULL);

	//ncextra = recv_fit_extra_constr_count(&fit);
	//if (ncextra)
	//	warning("Adding %zd %s to make parameters identifiable\n",
	//		ncextra, ncextra == 1 ? "constraint" : "constraints");

	//err = do_fit(&fit, NULL, NULL);


	//recv_fit_deinit(&fit);
terms_fail:
	design2_deinit(&d);
	design_deinit(&r);
	design_deinit(&s);
history_fail:
	history_deinit(&h);
properties_fail:
	error("Not implemented!");
	return NULL_USER_OBJECT;
}



static int get_properties(SEXP nsend, SEXP nrecv, SEXP loops,
			  struct properties *p)
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


	p->nsend = (size_t)xnsend;
	p->nrecv = (size_t)xnrecv;
	p->exclude_loops = !xloops;

	return 0;
}


static int get_history(SEXP time, SEXP sender, SEXP receiver, struct history *h)
{
	size_t nsend, nrecv;
	double *xtime;
	int *xsender, *xr;
	SEXP r;
	int i, n, j, m;
	struct message msg;
	size_t nto_max;
	size_t ito;


	/* extract dimensions */
	nsend = history_send_count(h);
	nrecv = history_recv_count(h);
	n = GET_LENGTH(time);


	/* setup message buffer */
	nto_max = nrecv;
	msg.to = (size_t *)R_alloc(nto_max, sizeof(size_t));
	msg.attr = 0;


	/* validate and extract 'time' */
	if (!IS_NUMERIC(time))
		DOMAIN_ERROR("'time' should be a numeric vector");

	xtime = NUMERIC_POINTER(time);


	/* validate and extract 'sender' */
	if (!IS_INTEGER(sender))
		DOMAIN_ERROR("'sender' should be an integer vector");
	if (GET_LENGTH(sender) != n)
		DOMAIN_ERROR("'time' and 'receiver' lengths differ");

	xsender = INTEGER_POINTER(sender);


	/* validate and extract 'receiver' */
	if (!IS_VECTOR(receiver))
		DOMAIN_ERROR("'receiver' should be a list");
	if (GET_LENGTH(receiver) != n)
		DOMAIN_ERROR("'time' and 'receiver' lengths differ");


	/* add all messages */
	for (i = 0; i < n; i++) {
		msg.time = xtime[i];
		msg.from = (size_t)(xsender[i] - 1);
		msg.attr = 0;

		/* validate and extract receiver[[i]] */
		r = VECTOR_ELT(receiver, i);
		if (!IS_INTEGER(r))
			DOMAIN_ERROR("each element of 'receiver' should be an integer vector");

		m = GET_LENGTH(r);
		xr = INTEGER_POINTER(r);

		for (j = 0; j < MIN(m, nto_max); j++) {
			msg.to[j] = (size_t)(xr[j] - 1);
		}
		msg.nto = (size_t)m;


		/* validate message */
		if (!R_FINITE(msg.time))
			DOMAIN_ERROR("'time' value is missing or infinite");
		if (msg.from >= nsend)
			DOMAIN_ERROR("'sender' value is out of range");
		if (msg.nto == 0)
			DOMAIN_ERROR("'receiver' value is empty");
		if (msg.nto > nrecv)
			DOMAIN_ERROR("'receiver' value contains duplicate elements");
		for (ito = 0; ito < msg.nto; ito++) {
			if (msg.to[ito] >= nrecv)
				DOMAIN_ERROR("'receiver' value is out of range");
		}


		/* add the message */
		history_set_time(h, msg.time);
		history_add(h, msg.from, msg.to, msg.nto, msg.attr);
	}

	return 0;
}


static int get_terms(SEXP factors, SEXP variables, SEXP order,
		     struct design *s, struct design *r, struct design2 *d,
		     struct terms *tm)
{
	return 0;
}

#if 0
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

#endif
