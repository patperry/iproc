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
	VARIABLE_TYPE_SEND_TRAIT,
	VARIABLE_TYPE_RECV_TRAIT,
	VARIABLE_TYPE_SEND_SPECIAL,
	VARIABLE_TYPE_RECV_SPECIAL,
	VARIABLE_TYPE_DYAD_SPECIAL
};


struct variable {
	enum variable_type type;
	const char *name;
	const char **names;
	size_t size;
	union {
		const struct var *send;
		const struct var *recv;
		const struct var2 *dyad;
	} var;
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


struct variables {
	struct variable *item;
	size_t count;
};



struct terms {
	struct variables variables;
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
static int get_terms(SEXP factors, SEXP term_labels, SEXP variables, SEXP order,
		     struct design *s, struct design *r, struct design2 *d,
		     struct terms *tm);
static int get_variables(SEXP variables,
			 struct design *s, struct design *r, struct design2 *d,
			 struct variables *v);
static int get_variable(SEXP variable, const char *name,
			struct design *s, struct design *r, struct design2 *d,
			struct variable *v);


static int do_fit(struct recv_fit *fit, const struct recv_params *params0, const double *duals0);


//static int get_ids(size_t *dst, SEXP src);



SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver,
		       SEXP factors, SEXP term_labels, SEXP variables, SEXP order,
		       SEXP nsend, SEXP nrecv, SEXP loops)
{
	struct properties p;
	struct history h;
	struct design s, r;
	struct design2 d;
	struct terms tm;

	const struct message *msgs;
	size_t nmsg;
	size_t ncextra;
	struct recv_fit fit;
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

	err = get_terms(factors, term_labels, variables, order, &s, &r, &d, &tm);
	if (err < 0)
		goto terms_fail;

	history_get_messages(&h, &msgs, &nmsg);

	recv_fit_init(&fit, &r, &d, p.exclude_loops, msgs, nmsg, NULL, NULL);

	ncextra = recv_fit_extra_constr_count(&fit);
	if (ncextra)
		warning("Adding %zd %s to make parameters identifiable\n",
			ncextra, ncextra == 1 ? "constraint" : "constraints");

	err = do_fit(&fit, NULL, NULL);

	recv_fit_deinit(&fit);
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
	if (!IS_INTEGER(nsend) || LENGTH(nsend) != 1)
		DOMAIN_ERROR("'nsend' should be a single integer");

	xnsend = INTEGER_VALUE(nsend);

	if (xnsend <= 0)
		DOMAIN_ERROR("'nsend' should be positive");


	/* validate and extract 'nrecv' */
	if (!IS_INTEGER(nrecv) || LENGTH(nrecv) != 1)
		DOMAIN_ERROR("'nrecv' should be a single integer");

	xnrecv = INTEGER_VALUE(nrecv);

	if (xnrecv <= 0)
		DOMAIN_ERROR("'nrecv' should be positive");



	/* validate and extract 'loops' */
	if (!IS_LOGICAL(loops) || LENGTH(loops) != 1)
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
	n = LENGTH(time);


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
	if (LENGTH(sender) != n)
		DOMAIN_ERROR("'time' and 'receiver' lengths differ");

	xsender = INTEGER_POINTER(sender);


	/* validate and extract 'receiver' */
	if (!IS_VECTOR(receiver))
		DOMAIN_ERROR("'receiver' should be a list");
	if (LENGTH(receiver) != n)
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

		m = LENGTH(r);
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


static int get_terms(SEXP factors, SEXP term_labels, SEXP variables, SEXP order,
		     struct design *s, struct design *r, struct design2 *d,
		     struct terms *tm)
{
	int i, j, m, n;
	int *xfactors;
	int nvariable, nterm;
	int err = 0;

	/* validate types */
	if (!IS_INTEGER(factors))
		DOMAIN_ERROR("'factors' should be integer");
	if (!isMatrix(factors))
		DOMAIN_ERROR("'factors' should be a matrix");
	if (!IS_CHARACTER(term_labels))
		DOMAIN_ERROR("'term.labels' should be a character vector");
	if (!IS_INTEGER(order))
		DOMAIN_ERROR("'order' should be a list");

	/* get dimensions */
	nterm = LENGTH(term_labels);
	nvariable = LENGTH(variables);

	/* validate factor matrix */
	if (isMatrix(factors)) {
		SEXP dim = GET_DIM(factors);
		int xf;

		m = INTEGER(dim)[0];
		n = INTEGER(dim)[1];
		xfactors = INTEGER_POINTER(factors);

		if ((size_t)m != nvariable)
			DOMAIN_ERROR("'factors' should have"
				     " rows equal to the number of variables");
		if ((size_t)n != nterm)
			DOMAIN_ERROR("'factors' should have"
				     " cols equal to the number of terms");

		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
				xf = xfactors[i * n + j];
				if (xf != 0 && xf != 1)
					DOMAIN_ERROR("invalid entry in 'factors' matrix");
			}
		}
	} else {
		m = 0;
		n = 0;
		xfactors = NULL;
	}

	err = get_variables(variables, s, r, d, &tm->variables);
	if (err < 0)
		goto out;

out:
	return err;
}


static int get_variables(SEXP variables,
			 struct design *s, struct design *r, struct design2 *d,
			 struct variables *v)
{
	SEXP names;
	int i, n;
	int err = 0;

	if (!IS_VECTOR(variables))
		DOMAIN_ERROR("'variables' should be a list");

	n = LENGTH(variables);

	names = GET_NAMES(variables);
	if (!IS_CHARACTER(names) && LENGTH(names) == n)
		DOMAIN_ERROR("'variables' must have a name for each component");

	v->count = (size_t)n;
	v->item = (void *)R_alloc(n, sizeof(*v->item));
	for (i = 0; i < n; i++) {
		SEXP var = VECTOR_ELT(variables, i);
		SEXP nm = STRING_ELT(names, i);

		err = get_variable(var, CHAR(nm), s, r, d, &v->item[i]);
		if (err < 0)
			goto out;
	}
out:
	return err;
}


static int get_variable(SEXP variable, const char *name,
			struct design *s, struct design *r, struct design2 *d,
			struct variable *v)
{
	size_t dims[VAR_RANK_MAX];
	size_t rank;

	double *xt, *x;
	int j, n, p;
	SEXP dim, cn, dn, nm;

	v->name = name;

	if (isMatrix(variable) || inherits(variable, "matrix")) {
		dim = GET_DIM(variable);
		n = INTEGER(dim)[0];
		p = INTEGER(dim)[1];
		xt = NUMERIC_POINTER(variable);
		x = (void *)R_alloc(n * p, sizeof(*x));

		if (p == 0)
			DOMAIN_ERROR("variable dimensions must be nonzero");

		if (n > 0)
			matrix_dtrans(n, p, xt, n, x, p);

		if (p == 1) {
			rank = 0;
		} else {
			rank = 1;
			dims[0] = (size_t)p;
		}

		if (inherits(variable, "send")) {
			if (n != design_count(s))
				DOMAIN_ERROR("variable dimension must match sender count");

			v->type = VARIABLE_TYPE_SEND_TRAIT;
			v->var.send = design_add_trait(s, NULL, x, dims, rank);
		} else {
			if (n != design_count(r))
				DOMAIN_ERROR("variable dimension must match receiver count");

			v->type = VARIABLE_TYPE_RECV_TRAIT;
			v->var.recv = design_add_trait(r, NULL, x, dims, rank);
		}

		dn = GET_DIMNAMES(variable);
		cn = GET_COLNAMES(dn);
		v->names = (void *)R_alloc(p, sizeof(*v->names));
		for (j = 0; j < p; j++) {
			nm = STRING_ELT(cn, j);
			v->names[j] = CHAR(nm);
		}

		v->size = (size_t)p;
	} else {
		if (inherits(variable, "send")) {
			v->type = VARIABLE_TYPE_SEND_SPECIAL;
		} else if (inherits(variable, "actor")) {
			v->type = VARIABLE_TYPE_RECV_SPECIAL;
		} else if (inherits(variable, "dyad")) {
			v->type = VARIABLE_TYPE_DYAD_SPECIAL;
		} else {
			DOMAIN_ERROR("unknown variable type");
		}

		v->size = 0;
	}
	return 0;
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


#if 0

static int get_ids(size_t *dst, SEXP src)
{
	int i, n;
	int *xsrc;
	int err = 0;

	xsrc = INTEGER_POINTER(src);
	n = LENGTH(src);

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
