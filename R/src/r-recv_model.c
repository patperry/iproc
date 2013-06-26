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


struct properties {
	/* dimensions */
	size_t nsend;
	size_t nrecv;

	/* properties */
	int exclude_loops;
	double skip;
};


enum variable_type {
	VARIABLE_TYPE_SEND_TRAIT,
	VARIABLE_TYPE_SEND_SPECIAL,
	VARIABLE_TYPE_RECV_TRAIT,
	VARIABLE_TYPE_RECV_SPECIAL,
	VARIABLE_TYPE_DYAD_TRAIT,
	VARIABLE_TYPE_DYAD_SPECIAL
};


struct variable {
	enum variable_type type;
	const char **names;
	union {
		const struct var *send;
		const struct var *recv;
		const struct var2 *dyad;
	} var;
};


struct variables {
	struct variable *item;
	size_t count;
};


struct terms {
	struct variable *item;
	size_t count;
};


struct terms_object {
	struct variables variables;
	struct terms terms;
};




static int get_properties(SEXP nsend, SEXP nrecv, SEXP loops, SEXP skip,
			  struct properties *p);
static int get_history(SEXP time, SEXP sender, SEXP receiver,
		       struct history *h);
static int get_terms_object(SEXP factors, SEXP term_labels, SEXP variables,
			    struct design *s, struct design *r, struct design2 *d,
			    struct terms_object *tm);
static int get_variables(SEXP variables,
			 struct design *s, struct design *r, struct design2 *d,
			 struct variables *v);
static int get_variable(SEXP variable, struct design *s, struct design *r, struct design2 *d,
			struct variable *v);
static int get_trait(SEXP variable, struct design *s, struct design *r, struct variable *v);
static int get_special(SEXP variable, struct design *s, struct design *r, struct design2 *d,
		       struct variable *v);
static int get_terms(SEXP factors, const struct variables *v,
		     struct design *s, struct design *r, struct design2 *d,
		     struct terms *tm);
static int get_term_factors(const int *factors, const struct variables *v, size_t j,
			    size_t *ind, size_t *n);
static int get_product(const struct variable *u, const struct variable *v,
		       struct variable *tm);


static int do_fit(struct recv_fit *fit, const struct recv_params *params0,
		  const double *duals0);





SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver,
		       SEXP factors, SEXP term_labels, SEXP variables,
		       SEXP nsend, SEXP nrecv, SEXP loops, SEXP skip)
{
	struct properties p;
	struct history h;
	struct design s, r;
	struct design2 d;
	struct terms_object tm;

	const struct message *msgs;
	size_t nmsg;
	size_t ncextra;
	double t0;
	struct recv_fit fit;
	int err = 0;

	err = get_properties(nsend, nrecv, loops, skip, &p);
	if (err < 0)
		goto properties_fail;

	history_init(&h, p.nsend, p.nrecv);
	err = get_history(time, sender, receiver, &h);
	if (err < 0)
		goto history_fail;

	design_init(&s, &h, p.nsend);
	design_init(&r, &h, p.nrecv);
	design2_init(&d, &h, p.nsend, p.nrecv);

	err = get_terms_object(factors, term_labels, variables, &s, &r, &d, &tm);
	if (err < 0)
		goto terms_fail;

	history_get_messages(&h, &msgs, &nmsg);

	t0 = nmsg ? msgs[0].time : -INFINITY;
	while (nmsg && msgs->time < t0 + p.skip) {
		nmsg--;
		msgs++;
	}

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



static int get_properties(SEXP nsend, SEXP nrecv, SEXP loops, SEXP skip,
			  struct properties *p)
{
	int xnsend, xnrecv, xloops;
	double xskip;

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


	/* validate and extract 'skip' */
	if (!IS_NUMERIC(skip) || LENGTH(skip) != 1)
		DOMAIN_ERROR("'skip' should be a single numeric value");

	xskip = NUMERIC_VALUE(skip);
	if (!(xskip >= 0))
		DOMAIN_ERROR("'skip' should be non-negative");


	p->nsend = (size_t)xnsend;
	p->nrecv = (size_t)xnrecv;
	p->exclude_loops = !xloops;
	p->skip = xskip;

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


static int get_terms_object(SEXP factors, SEXP term_labels, SEXP variables,
			    struct design *s, struct design *r, struct design2 *d,
			    struct terms_object *tm)
{
	int err = 0;

	/* validate term.labels */
	if (!IS_CHARACTER(term_labels))
		DOMAIN_ERROR("'term.labels' should be a character vector");


	err = get_variables(variables, s, r, d, &tm->variables);
	if (err < 0)
		goto out;

	err = get_terms(factors, &tm->variables, s, r, d, &tm->terms);
	if (err < 0)
		goto out;

out:
	return err;
}


static int get_variables(SEXP variables,
			 struct design *s, struct design *r, struct design2 *d,
			 struct variables *v)
{
	SEXP var;
	int i, n;
	int err = 0;

	if (!IS_VECTOR(variables))
		DOMAIN_ERROR("'variables' should be a list");

	n = LENGTH(variables);

	v->count = (size_t)n;
	v->item = (void *)R_alloc(n, sizeof(*v->item));
	for (i = 0; i < n; i++) {
		var = VECTOR_ELT(variables, i);

		err = get_variable(var, s, r, d, &v->item[i]);
		if (err < 0)
			goto out;
	}
out:
	return err;
}


static int get_variable(SEXP variable, struct design *s, struct design *r, struct design2 *d,
			struct variable *v)
{
	int err = 0;

	if (isMatrix(variable) || inherits(variable, "matrix")) {
		err = get_trait(variable, s, r, v);
	} else if (inherits(variable, "special")) {
		err = get_special(variable, s, r, d, v);
	} else {
		DOMAIN_ERROR("unknown variable type");
	}

	return err;
}

		
static int get_trait(SEXP variable, struct design *s, struct design *r, struct variable *v)
{
	SEXP dim, cn, dn, nm;
	int j, n, p;
	double *xt, *x;
	size_t dims[VAR_RANK_MAX];
	size_t rank;

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

	return 0;
}


static int get_special(SEXP variable, struct design *s, struct design *r, struct design2 *d,
		       struct variable *v)
{
	SEXP i1, i2, labels, l;
	double xi, *xi1, *xi2;
	size_t ni1, ni2, k, i, j, size, rank;
	struct design *a;
	const struct var **var;
	const struct var2 **var2;

	if (!(inherits(variable, "interval")
		|| inherits(variable, "intervals")
		|| inherits(variable, "intervals2")))
		DOMAIN_ERROR("unknown variable type");


	/* extract intervals */
	i1 = VECTOR_ELT(variable, 0);
	xi1 = NUMERIC_POINTER(i1);
	ni1 = LENGTH(i1);

	if (inherits(variable, "interval")) {
		xi = *xi1;
		size = 1;
		rank = 0;
	} else if (inherits(variable, "intervals")) {
		size = ni1;
		rank = 1;
	} else {
		/* intervals2 */
		i2 = VECTOR_ELT(variable, 1);
		xi2 = NUMERIC_POINTER(i2);
		ni2 = LENGTH(i2);
		size = ni1 * ni2;
		rank = 2;
	}


	/* set variable names */
	labels = getAttrib(variable, install("labels"));
	v->names = (void *)R_alloc(size, sizeof(*v->names));
	for (k = 0; k < size; k++) {
		l = STRING_ELT(labels, k);

		if (rank == 2) {
			i = k % ni1;
			j = k / ni1;
			v->names[i * ni2 + j] = CHAR(l);
		} else {
			v->names[k] = CHAR(l);
		}
	}


	/* set type */
	if (inherits(variable, "send")) {
		a = s;
		v->type = VARIABLE_TYPE_SEND_SPECIAL;
		var = &v->var.send;
	} else if (inherits(variable, "actor")) {
		a = r;
		v->type = VARIABLE_TYPE_RECV_SPECIAL;
		var = &v->var.recv;
	} else if (inherits(variable, "dyad")) {
		v->type = VARIABLE_TYPE_DYAD_SPECIAL;
		var2 = &v->var.dyad;
	} else {
		DOMAIN_ERROR("unknown variable type");
	}


	/* add variable to design */
	if (inherits(variable, "irecvtot")) {
		*var = design_add_tvar(a, NULL, VAR_IRECVTOT, xi);
	} else if (inherits(variable, "isendtot")) {
		*var = design_add_tvar(a, NULL, VAR_ISENDTOT, xi);
	} else if (inherits(variable, "nrecvtot")) {
		*var = design_add_tvar(a, NULL, VAR_NRECVTOT, xi1, ni1);
	} else if (inherits(variable, "nsendtot")) {
		*var = design_add_tvar(a, NULL, VAR_NSENDTOT, xi1, ni1);
	} else if (inherits(variable, "irecv")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_IRECV, xi);
	} else if (inherits(variable, "isend")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_ISEND, xi);
	} else if (inherits(variable, "nrecv")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_NRECV, xi1, ni1);
	} else if (inherits(variable, "nsend")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_NSEND, xi1, ni1);
	} else if (inherits(variable, "nrecv2")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_NRECV2, xi1, ni1, xi2, ni2);
	} else if (inherits(variable, "nsend2")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_NSEND2, xi1, ni1, xi2, ni2);
	} else if (inherits(variable, "nsib")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_NSIB, xi1, ni1, xi2, ni2);
	} else if (inherits(variable, "ncosib")) {
		*var2 = design2_add_tvar(d, NULL, VAR2_NCOSIB, xi1, ni1, xi2, ni2);
	} else {
		DOMAIN_ERROR("unknown variable type");
	}

	return 0;
}


static int get_terms(SEXP factors, const struct variables *v,
		     struct design *s, struct design *r, struct design2 *d,
		     struct terms *tm)
{
	SEXP dim;
	int j, m, n, *xfactors;
	size_t *ind, order;
	int err = 0;

	if (!IS_INTEGER(factors))
		DOMAIN_ERROR("'factors' should be integer");
	if (!isMatrix(factors))
		DOMAIN_ERROR("'factors' should be a matrix");


	/* validate factor matrix */
	if (isMatrix(factors)) {
		dim = GET_DIM(factors);
		m = INTEGER(dim)[0];
		n = INTEGER(dim)[1];
		xfactors = INTEGER_POINTER(factors);

		if ((size_t)m != v->count)
			DOMAIN_ERROR("'factors' should have"
				     " rows equal to the number of variables");
	} else {
		m = 0;
		n = 0;
		xfactors = NULL;
	}

	ind = (void *)R_alloc(v->count, sizeof(*ind));


	tm->item = (void *)R_alloc(n, sizeof(*tm->item));
	tm->count = n;

	for (j = 0; j < n; j++) {
		err = get_term_factors(xfactors, v, j, ind, &order);
		if (err < 0)
			goto out;

		if (order == 0) {
			DOMAIN_ERROR("intercepts are not allowed");
		} else if (order == 1) {
			memcpy(&tm->item[j], &v->item[ind[0]],
			       sizeof(struct variable));
		} else if (order == 2) {
			err = get_product(v->item + ind[0], v->item + ind[1],
					  tm->item + j);
			if (err < 0)
				goto out;
		} else {
			DOMAIN_ERROR("interactions of order >= 3 are not implemented yet");
		}

	}

out:
	return err;
}


static int get_term_factors(const int *factors, const struct variables *v, size_t j,
			    size_t *ind, size_t *n)
{
	size_t i, m, order;
	int f;

	m = v->count;
	order = 0;

	for (i = 0; i < m; i++) {
		f = factors[i + j * m];
		if (f == 1) {
			ind[order++] = i;
		} else if (f != 0) {
			DOMAIN_ERROR("invalid entry in 'factors' matrix");
		}
	}

	*n = order;
	return 0;
}


static int get_product(const struct variable *u, const struct variable *v,
		       struct variable *tm)
{
	DOMAIN_ERROR("interactions are not implemented yet");
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

