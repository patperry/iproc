#include "port.h"

#include <assert.h>
#include <errno.h>
#include <stddef.h>
#include <stdio.h>
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

enum variable_design {
	VARIABLE_DESIGN_SEND,
	VARIABLE_DESIGN_RECV,
	VARIABLE_DESIGN_DYAD
};

enum variable_type {
	VARIABLE_TYPE_TRAIT,
	VARIABLE_TYPE_SPECIAL
};


struct variable {
	enum variable_design design;
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
static int get_terms_object(SEXP factors, SEXP variables, struct design *s,
		            struct design *r, struct design2 *d,
			    struct terms_object *tm);
static int get_variables(SEXP variables, struct design *s, struct design *r,
		         struct design2 *d, struct variables *v);
static int get_variable(SEXP variable, struct design *s, struct design *r,
			struct design2 *d, struct variable *v);
static int get_trait(SEXP variable, struct design *s, struct design *r,
		     struct variable *v);
static int get_special(SEXP variable, struct design *s, struct design *r,
		       struct design2 *d, struct variable *v);
static int get_terms(SEXP factors, const struct variables *v, struct design *s,
		     struct design *r, struct design2 *d, struct terms *tm);
static int get_term_factors(const char *name, const int *factors,
			    const struct variables *v, size_t j, size_t *ind,
			    size_t *n);
static int get_product(const char *name, const struct variable *u,
		       const struct variable *v, struct design *s,
		       struct design *r, struct design2 *d,
		       struct variable *tm);

static int do_fit(struct recv_fit *fit, const struct recv_params *params0,
		  const double *duals0, size_t *iter);


static SEXP alloc_term_labels(const struct terms *tm);
static size_t terms_size(const struct terms *tm);
static size_t variable_size(const struct variable *v);
static const struct var_meta *variable_meta(const struct variable *v);

static void get_variable_names(const struct variable *v, SEXP dst, size_t off);
static void cindex_to_coord(const size_t *dims, size_t rank, size_t i,
			    size_t *coord);
static size_t coord_to_findex(const size_t *dims, size_t rank,
			      const size_t *coord);

static SEXP alloc_terms_permute(const struct terms *tm,
				const struct design *r,
				const struct design2 *d);

static SEXP alloc_vector_copy(const double *x, size_t n);
static SEXP alloc_matrix_copy(const double *x, size_t m, size_t n);
static SEXP alloc_score(const struct recv_loglik *ll);
static SEXP alloc_imat(const struct recv_fit *fit);


SEXP Riproc_recv_model(SEXP time, SEXP sender, SEXP receiver, SEXP factors,
		       SEXP variables, SEXP nsend, SEXP nrecv, SEXP loops,
		       SEXP skip)
{
	struct properties p;
	struct history h;
	struct design s, r;
	struct design2 d;
	struct terms_object tm;
	const struct message *msgs;
	size_t nmsg;
	size_t iter, dim, nc, ntot, rank;
	size_t ncextra;
	double t0;
	struct recv_fit fit;
	const struct recv_loglik *ll;
	const struct recv_model *model;
	const struct constr *constr;
	const double *beta, *nu;
	double dev;
	SEXP names;
	SEXP ret = NULL_USER_OBJECT;
	int k;
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

	err = get_terms_object(factors, variables, &s, &r, &d, &tm);
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
		warning("adding %zd %s to make parameters identifiable\n",
			ncextra, ncextra == 1 ? "constraint" : "constraints");

	err = do_fit(&fit, NULL, NULL, &iter);
	if (err < 0)
		goto fit_fail;

	constr = recv_fit_constr(&fit);
	nc = constr_count(constr);
	ll = recv_fit_loglik(&fit);
	ntot = recv_loglik_count(ll);
	model = recv_loglik_model(ll);
	beta = (recv_model_params(model))->recv.traits; /* HACK */
	nu = recv_fit_duals(&fit);
	dim = recv_model_dim(model);
	rank = dim - nc;
	dev = recv_loglik_dev(ll);

	PROTECT(ret = NEW_LIST(16));
	PROTECT(names = NEW_CHARACTER(16));
	SET_NAMES(ret, names);
	k = 0;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("coefficients"));
	SET_VECTOR_ELT(ret, k, alloc_vector_copy(beta, dim));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("duals"));
	SET_VECTOR_ELT(ret, k, alloc_vector_copy(nu, nc));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("constraints"));
	SET_VECTOR_ELT(ret, k, alloc_matrix_copy(constr_all_wts(constr), nc, dim));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("constraint.values"));
	SET_VECTOR_ELT(ret, k, alloc_vector_copy(constr_all_vals(constr), nc));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("rank"));
	SET_VECTOR_ELT(ret, k, ScalarInteger((int)rank));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("deviance"));
	SET_VECTOR_ELT(ret, k, ScalarReal(dev));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("aic"));
	SET_VECTOR_ELT(ret, k, ScalarReal(dev + (double) 2 * rank));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("null.deviance"));
	SET_VECTOR_ELT(ret, k, ScalarReal(recv_fit_dev0(&fit)));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("iter"));
	SET_VECTOR_ELT(ret, k, ScalarInteger((int)iter));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("df.residual"));
	SET_VECTOR_ELT(ret, k, ScalarReal((double)(ntot - rank)));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("df.null"));
	SET_VECTOR_ELT(ret, k, ScalarReal((double)ntot));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("converged"));
	SET_VECTOR_ELT(ret, k, ScalarLogical(TRUE));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("names"));
	SET_VECTOR_ELT(ret, k, alloc_term_labels(&tm.terms));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("perm"));
	SET_VECTOR_ELT(ret, k, alloc_terms_permute(&tm.terms, &r, &d));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("score"));
	SET_VECTOR_ELT(ret, k, alloc_score(ll));
	k++;

	SET_STRING_ELT(names, k, COPY_TO_USER_STRING("imat"));
	SET_VECTOR_ELT(ret, k, alloc_imat(&fit));
	k++;

	UNPROTECT(2);


fit_fail:
	recv_fit_deinit(&fit);
terms_fail:
	design2_deinit(&d);
	design_deinit(&r);
	design_deinit(&s);
history_fail:
	history_deinit(&h);
properties_fail:
	return ret;
}


static SEXP alloc_vector_copy(const double *x, size_t n)
{
	SEXP dst = NEW_NUMERIC(n);
	memcpy(NUMERIC_POINTER(dst), x, n * sizeof(double));
	return dst;
}

static SEXP alloc_matrix_copy(const double *x, size_t m, size_t n)
{
	SEXP dst = allocMatrix(REALSXP, m, n);
	double *xdst = NUMERIC_POINTER(dst);

	if (m > 0 && n > 0)
		matrix_dtrans(n, m, x, n, xdst, m);

	return dst;
}

static SEXP alloc_score(const struct recv_loglik *ll)
{
	const struct recv_model *m = recv_loglik_model(ll);
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	struct recv_params p;
	size_t dim = recv_model_dim(m);
	SEXP score = NEW_NUMERIC(dim);
	double *xscore = NUMERIC_POINTER(score);

	p.recv.traits = xscore;
	p.recv.tvars = p.recv.traits + design_trait_dim(r);
	p.dyad.traits = p.recv.tvars + design_tvar_dim(r);
	p.dyad.tvars = p.dyad.traits + design2_trait_dim(d);

	memset(xscore, 0, dim * sizeof(double));
	recv_loglik_axpy_score(1.0, ll, &p);

	return score;
}


static SEXP alloc_imat(const struct recv_fit *fit)
{
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_model *m = recv_loglik_model(ll);
	size_t dim = recv_model_dim(m);
	const double *imatp;
	SEXP imat = allocMatrix(REALSXP, dim, dim);
	double *ximat = NUMERIC_POINTER(imat);
	enum blas_uplo uplo, f77uplo;
	size_t i, j;

	recv_fit_get_imat(fit, &imatp, &uplo);
	f77uplo = (uplo == BLAS_UPPER) ? BLAS_LOWER : BLAS_UPPER;

	memset(ximat, 0, dim * dim * sizeof(double));

	if (dim > 0) {
		packed_dsctr(f77uplo, dim, imatp, ximat, dim);
	}

	if (uplo == BLAS_UPPER) {
		for (j = 0; j < dim; j++) {
			for (i = 0; i < j; i++) {
				ximat[i + j * dim] = ximat[j + i * dim];
			}
		}
	} else {
		for (j = 0; j < dim; j++) {
			for (i = j + 1; i < dim; i++) {
				ximat[i + j * dim] = ximat[j + i * dim];
			}
		}
	}

	return imat;
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


static int get_terms_object(SEXP factors, SEXP variables, struct design *s,
			    struct design *r, struct design2 *d,
			    struct terms_object *tm)
{
	int err = 0;

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

		v->design = VARIABLE_DESIGN_SEND;
		v->var.send = design_add_trait(s, NULL, x, dims, rank);
	} else {
		if (n != design_count(r))
			DOMAIN_ERROR("variable dimension must match receiver count");

		v->design = VARIABLE_DESIGN_RECV;
		v->var.recv = design_add_trait(r, NULL, x, dims, rank);
	}
	v->type = VARIABLE_TYPE_TRAIT;

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
		v->design = VARIABLE_DESIGN_SEND;
		var = &v->var.send;
	} else if (inherits(variable, "actor")) {
		a = r;
		v->design = VARIABLE_DESIGN_RECV;
		var = &v->var.recv;
	} else if (inherits(variable, "dyad")) {
		v->design = VARIABLE_DESIGN_DYAD;
		var2 = &v->var.dyad;
	} else {
		DOMAIN_ERROR("unknown variable type");
	}
	v->type = VARIABLE_TYPE_SPECIAL;


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
	SEXP dim, dimnames, colnames;
	int j, m, n, *xfactors;
	size_t *ind, order;
	int err = 0;
	const char *name;

	if (!IS_INTEGER(factors))
		DOMAIN_ERROR("'factors' should be integer");
	if (!(isMatrix(factors) || LENGTH(factors) == 0))
		DOMAIN_ERROR("'factors' should be a matrix or an empty integer vector");


	/* validate factor matrix */
	if (isMatrix(factors)) {
		dim = GET_DIM(factors);
		m = INTEGER(dim)[0];
		n = INTEGER(dim)[1];
		xfactors = INTEGER_POINTER(factors);

		if ((size_t)m != v->count)
			DOMAIN_ERROR("'factors' should have"
				     " rows equal to the number of variables");

		dimnames = GET_DIMNAMES(factors);
		colnames = GET_COLNAMES(dimnames);
	} else {
		m = 0;
		n = 0;
		xfactors = NULL;
	}

	ind = (void *)R_alloc(v->count, sizeof(*ind));


	tm->item = (void *)R_alloc(n, sizeof(*tm->item));
	tm->count = n;

	for (j = 0; j < n; j++) {
		name = CHAR(STRING_ELT(colnames, j));
		err = get_term_factors(name, xfactors, v, j, ind, &order);
		if (err < 0)
			goto out;

		if (order == 0) {
			DOMAIN_ERROR("intercepts are not allowed");
		} else if (order == 1) {
			memcpy(&tm->item[j], &v->item[ind[0]],
			       sizeof(struct variable));
		} else if (order == 2) {
			err = get_product(name,
					  v->item + ind[0], v->item + ind[1],
					  s, r, d, tm->item + j);
			if (err < 0)
				goto out;
		} else {
			DOMAIN_ERROR("interactions of order >= 3 are not implemented yet");
		}

	}

out:
	return err;
}


static int get_term_factors(const char *name,
		            const int *factors, const struct variables *v, size_t j,
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
			error("invalid entry in 'factors' matrix for term '%s'", name);
			return -EDOM;
		}
	}

	*n = order;
	return 0;
}


static int get_product(const char *name, const struct variable *u, const struct variable *v,
		       struct design *s, struct design *r, struct design2 *d,
		       struct variable *tm)
{
	int err = -EDOM;
	size_t k, l, nu, nv, n, len, lenu, lenv;
	const char *stru, *strv;
	char *str;

	if (u->design == VARIABLE_DESIGN_SEND) {
		nu = u->var.send->meta.size;

		if (v->design == VARIABLE_DESIGN_RECV) {
			nv = v->var.recv->meta.size;
			tm->design = VARIABLE_DESIGN_DYAD;

			if (u->type == VARIABLE_TYPE_TRAIT
					&& v->type == VARIABLE_TYPE_TRAIT) {
				tm->type = VARIABLE_TYPE_TRAIT;
				tm->var.dyad = design2_add_kron(d, NULL,
								u->var.send,
								v->var.recv);
				err = 0;
			}
			/* TODO s(trait):r(special) */
			/* TODO s(special):r(trait) */
			/* TODO s(special):r(special) */
		}
		/* TODO s(.):d(.) */
	} else if (u->design == VARIABLE_DESIGN_RECV) {
		nu = u->var.recv->meta.size;

		if (v->design == VARIABLE_DESIGN_RECV) {
			nv = v->var.recv->meta.size;
			tm->design = VARIABLE_DESIGN_RECV;

			if (u->type == VARIABLE_TYPE_TRAIT
					&& v->type == VARIABLE_TYPE_TRAIT) {
				tm->type = VARIABLE_TYPE_TRAIT;
			} else {
				tm->type = VARIABLE_TYPE_SPECIAL;
			}

			tm->var.recv = design_add_prod(r, NULL, u->var.recv,
						       v->var.recv);
			err = 0;
		}
		/* TODO r(.):d(.) */

	} else if (u->design == VARIABLE_DESIGN_DYAD) {
		nu = u->var.dyad->meta.size;

		if (v->design == VARIABLE_DESIGN_DYAD) {
			nv = v->var.dyad->meta.size;
			tm->design = VARIABLE_DESIGN_DYAD;

			if (u->type != VARIABLE_TYPE_TRAIT
					|| v->type != VARIABLE_TYPE_TRAIT) {
				tm->type = VARIABLE_TYPE_SPECIAL;

				tm->var.dyad = design2_add_prod(d, NULL,
								u->var.dyad,
								v->var.dyad);
				err = 0;
			}
			/* TODO d(trait):d(trait) */
		}
	}

	if (!err) {
		n = nu * nv;
		tm->names = (void *)R_alloc(n, sizeof(*tm->names));

		for (k = 0; k < nu; k++) {
			stru = u->names[k];
			lenu = strlen(stru);

			for (l = 0; l < nv; l++) {
				strv = v->names[l];
				lenv = strlen(strv);

				len = lenu + lenv + 2;
				str = (void *)R_alloc(len, sizeof(*str));


				sprintf(str, "%s:%s", stru, strv);
				tm->names[k * nv + l] = str;
			}
		}

		/*Rprintf("`%s'\n", name);
		for (k = 0; k < n; k++) {
			Rprintf("\t`%s'\n", tm->names[k]);
		}
		*/
	}


	if (err < 0)
		error("interaction for terms of type '%s' are not implemented yet", name);

	return err;
}


static int do_fit(struct recv_fit *fit, const struct recv_params *params0,
		  const double *duals0, size_t *iter)
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

	*iter = it;

	return err;
}


static SEXP alloc_term_labels(const struct terms *tm)
{
	SEXP labels;
	size_t ind, len, i, m, n;
	const struct variable *v;

	len = terms_size(tm);
	m = tm->count;

	PROTECT(labels = NEW_STRING(len));
	ind = 0;
	for (i = 0; i < m; i++) {
		v = &tm->item[i];
		n = variable_size(v);
		get_variable_names(v, labels, ind);
		ind += n;
	}
	UNPROTECT(1);

	return labels;
}




static size_t terms_size(const struct terms *tm)
{
	size_t size, i, n;

	size = 0;
	n = tm->count;
	for (i = 0; i < n; i++) {
		size += variable_size(&tm->item[i]);
	}

	return size;
}


static SEXP alloc_terms_permute(const struct terms *tm,
				const struct design *r,
				const struct design2 *d)
{
	size_t i, j, jf, ind, m, n, len;
	const struct variable *v;
	const struct var_meta *meta;
	size_t coord[VAR_RANK_MAX];
	size_t dimr0, dimr, dimd0;
	size_t off;
	SEXP perm;
	int *xperm;

	len = terms_size(tm);
	PROTECT(perm = NEW_INTEGER(len));
	xperm = INTEGER_POINTER(perm);

	dimr0 = design_trait_dim(r);
	dimr = design_dim(r);
	dimd0 = design2_trait_dim(d);

	m = tm->count;

	ind = 0;
	for (i = 0; i < m; i++) {
		v = &tm->item[i];
		meta = variable_meta(v);

		if (v->design == VARIABLE_DESIGN_RECV) {
			if (v->type == VARIABLE_TYPE_TRAIT) {
				off = 0;
			} else {
				off = dimr0;
			}
			off += v->var.recv->index;
		} else {
			if (v->type == VARIABLE_TYPE_TRAIT) {
				off = dimr;
			} else {
				off = dimr + dimd0;
			}
			off += v->var.dyad->index;
		}

		n = meta->size;
		for (j = 0; j < n; j++) {
			cindex_to_coord(meta->dims, meta->rank, j, coord);
			jf = coord_to_findex(meta->dims, meta->rank, coord);
			xperm[ind + jf] = (int)(off + j + 1);
		}
		ind += n;
	}

	UNPROTECT(1);
	return perm;
}


static size_t variable_size(const struct variable *v)
{
	const struct var_meta *meta = variable_meta(v);
	return meta->size;
}


static const struct var_meta *variable_meta(const struct variable *v)
{
	switch (v->design) {
	case VARIABLE_DESIGN_SEND:
		return &v->var.send->meta;
	case VARIABLE_DESIGN_RECV:
		return &v->var.recv->meta;
	case VARIABLE_DESIGN_DYAD:
		return &v->var.dyad->meta;
	}

	assert(0);
	return NULL;
}


static void get_variable_names(const struct variable *v, SEXP dst, size_t off)
{
	const struct var_meta *meta = variable_meta(v);
	const char **names = v->names;
	size_t i, j, n = meta->size;
	size_t coord[VAR_RANK_MAX];

	for (i = 0; i < n; i++) {
		cindex_to_coord(meta->dims, meta->rank, i, coord);
		j = coord_to_findex(meta->dims, meta->rank, coord);

		SET_STRING_ELT(dst, off + j, COPY_TO_USER_STRING(names[i]));
	}
}


static void cindex_to_coord(const size_t *dims, size_t rank, size_t i,
			    size_t *coord)
{
	/*
	 * array1: i1
	 * array2: i1 *  n2            + i2
	 * array3: i1 * (n2 * n3)      + i2 *  n3       + i3
	 * array4: i1 * (n2 * n3 * n4) + i2 * (n3 * n4) + i3 * n4 + i4
	 *
	 **/

	if (rank == 1) {
		coord[0] = i;
	} else if (rank > 1) {
		coord[rank - 1] = i % dims[rank - 1];
		cindex_to_coord(dims, rank - 1, i / dims[rank - 1], coord);
	}
}


static size_t coord_to_findex(const size_t *dims, size_t rank,
			      const size_t *coord)
{
	size_t r, i, stride;

	stride = 1;
	i = 0;
	for (r = 0; r < rank; r++) {
		i += coord[r] * stride;
		stride *= dims[r];
	}

	return i;
}

