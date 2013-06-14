#include "port.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "coreutil.h"
#include "xalloc.h"

#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "dSFMT/dSFMT.h"

#include "enron/actors.h"
#include "enron/messages.h"

#include "recv_boot.h"
#include "recv_fit.h"
#include "recv_resid.h"


static int setup_history(struct history *h)
{
	size_t nsend = ENRON_ACTOR_COUNT;
	size_t nrecv = nsend;
	size_t maxrecip = 5; /* exclude messages with more than 5 recipients */

	double *time;
	size_t *from;
	size_t **to;
	size_t *nto;
	ptrdiff_t *attr;
	size_t i, nmsg;
	int err;

	err = enron_messages_init(maxrecip, &time, &from, &to, &nto, &attr, &nmsg);
	if (err)
		return err;

	history_init(h, nsend, nrecv);
	for (i = 0; i < nmsg; i++) {
		history_set_time(h, time[i]);
		history_add(h, from[i], to[i], nto[i], attr[i]);
		free(to[i]);
	}

	free(attr);
	free(nto);
	free(to);
	free(from);
	free(time);

	return err;
}


static int setup_designs(struct design *r, struct design2 *d, struct history *h)
{
	size_t nsend = history_nsend(h);
	size_t nrecv = history_nrecv(h);

	double intvls[] = {
		/* 450, */	/*   7.5m */
		/* 900, */	/*  15.0m */
		1800,		/*  30.0m */
		/* 3600, */	/*   1.0h */
		7200,		/*   2.0h */
		/* 14400, */	/*   4.0h */
		28800,		/*   8.0h */
		/* 57600, */	/*  16.0h */
		115200,		/*  32.0h */
		/* 230400, */	/*   2.6d */
		460800,		/*   5.3d */
		/* 921600, */	/*  10.6d */
		1843200,	/*  21.3d */
		/* 3686400, */	/*  42.6d */
		/* 7372800, */	/*  85.3d */
		/* 14745600, */	/* 170.6d */
		INFINITY	/*   Inf  */
	};
	size_t nintvl = sizeof(intvls) / sizeof(intvls[0]);
	double window = INFINITY;

	size_t terms = 1; /* first-order interactions only */
	double *trait_x;
	const char * const *trait_names;
	size_t trait_dim;
	int err;

	char namebuf[1024];
	size_t k, l;


	err = enron_actors_init(terms, &trait_x, &trait_names, &trait_dim);
	if (err)
		return err;

	design_init(r, h, nsend);
	design2_init(d, h, nsend, nrecv);


	/* recv design */
	design_add_traits(r, trait_names, trait_x, trait_dim);
	//design_add_prod(r, "Leg:Jun", design_var(r, "Leg"), design_var(r, "Jun"));
	//design_add_prod(r, "Trad:Jun", design_var(r, "Trad"), design_var(r, "Jun"));
	//design_add_prod(r, "Leg:Fem", design_var(r, "Leg"), design_var(r, "Fem"));
	//design_add_prod(r, "Trad:Fem", design_var(r, "Trad"), design_var(r, "Fem"));
	//design_add_prod(r, "Jun:Fem", design_var(r, "Jun"), design_var(r, "Fem"));

	/* third order */
	//design_add_prod(r, "Leg:Jun:Fem", design_var(r, "Leg:Jun"), design_var(r, "Fem"));
	//design_add_prod(r, "Trad:Jun:Fem", design_var(r, "Trad:Jun"), design_var(r, "Fem"));

	//design_add_tvar(r, "ISendTot", VAR_ISENDTOT, window);
	//design_add_tvar(r, "IRecvTot", VAR_IRECVTOT, window);
	//design_add_prod(r, "ISendTot:IRecvTot", design_var(r, "ISendTot"), design_var(r, "IRecvTot"));

	//design_add_tvar(r, "NRecvTot", VAR_NRECVTOT, intvls, nintvl);
	//design_add_tvar(r, "NSendTot", VAR_NSENDTOT, intvls, nintvl);
	//design_add_prod(r, "NRecvTot:Fem", design_var(r, "NRecvTot"), design_var(r, "Fem"));


	/* dyad design */
	//design2_add_tvar(d, "ISend", VAR2_ISEND, window);
	//design2_add_tvar(d, "IRecv", VAR2_IRECV, window);
	//design2_add_prod(d, "ISend:IRecv", design2_var(d, "ISend"), design2_var(d, "IRecv"));

	//design2_add_tvar(d, "NSend", VAR2_NSEND, intvls, nintvl);
	//design2_add_tvar(d, "NRecv", VAR2_NRECV, intvls, nintvl);

	//design2_add_tvar(d, "NSend2", VAR2_NSEND2, intvls, nintvl, intvls, nintvl);
	//design2_add_tvar(d, "NRecv2", VAR2_NRECV2, intvls, nintvl, intvls, nintvl);
	//design2_add_tvar(d, "NSib", VAR2_NSIB, intvls, nintvl, intvls, nintvl);
	//design2_add_tvar(d, "NCosib", VAR2_NCOSIB, intvls, nintvl, intvls, nintvl);

	/*for (k = 0; k < design_trait_count(r); k++) {
		const struct var *u = design_trait_item(r, k);
		for (l = 0; l < design_trait_count(r); l++) {
			const struct var *v = design_trait_item(r, l);
			sprintf(namebuf, "%s*%s", u->meta.name, v->meta.name);
			design2_add_kron(d, namebuf, u, v);
		}
	}
	 */

	free(trait_x);

	return err;
}


static void setup_constr(struct constr *c, const struct design *r, const struct design2 *d)
{
	size_t dim = design_dim(r) + design2_dim(d);
	constr_init(c, dim);
}


/*

#define INTERVALS		"intervals"
#define VARIATE_NAMES		"variate_names"
#define COHORT_NAMES		"cohort_names"
#define COEFFICIENTS		"coefficients"
#define COUNT			"count"
#define SCORE			"score"
#define INFORMATION		"information"
#define RANK			"rank"
#define CONSTRAINTS		"constraints"
#define CONSTRAINT_VALUES	"constraint_values"
#define DUALS			"duals"
#define DEVIANCE		"deviance"
#define NULL_DEVIANCE		"null_deviance"
#define DF_RESIDUAL		"df_residual"
#define DF_NULL			"df_null"
#define FITTED_COUNTS		"fitted_counts"
#define OBSERVED_COUNTS		"observed_counts"



static const char **alloc_varnames(const struct frame *f, size_t *pn)
{
	const struct var_name_fmt fmt = VAR_NAME_FMT0;
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t index = 0, dim = recv_coefs_dim(f);
	size_t k, n, i;

	const char **res = xmalloc(dim * sizeof(*res));

	n = design_trait_count(r);
	for (k = 0; k < n; k++) {
		const struct var_meta *v = &(design_trait_var(r, k))->meta;
		for (i = 0; i < v->size; i++) {
			res[index++] = alloc_var_name(&fmt, v, i);
		}
	}

	n = design_tvar_count(r);
	for (k = 0; k < n; k++) {
		const struct var_meta *v = &(design_tvar_var(r, k))->meta;
		for (i = 0; i < v->size; i++) {
			res[index++] = alloc_var_name(&fmt, v, i);
		}
	}

	n = design2_trait_count(d);
	for (k = 0; k < n; k++) {
		const struct var_meta *v = &(design2_trait_var(d, k))->meta;
		for (i = 0; i < v->size; i++) {
			res[index++] = alloc_var_name(&fmt, v, i);
		}
	}

	n = design2_tvar_count(d);
	for (k = 0; k < n; k++) {
		const struct var_meta *v = &(design2_tvar_var(d, k))->meta;
		for (i = 0; i < v->size; i++) {
			res[index++] = alloc_var_name(&fmt, v, i);
		}
	}
	assert(index == dim);

	*pn = dim;

	return res;
}

static void free_varnames(const char **varnames, size_t nvarnames)
{
	size_t i;

	for (i = 0; i < nvarnames; i++) {
		free((char *)(varnames[i]));
	}
	free(varnames);
}


static herr_t output_variate_names(hid_t loc_id, const struct frame *f)
{
	hid_t dataset_id, dataspace_id;
	herr_t status;

	hsize_t dims[1];
	size_t nvarnames;
	const char **varnames = alloc_varnames(f, &nvarnames);

	hid_t type;
	type = H5Tcopy(H5T_C_S1);
	status = H5Tset_size(type, H5T_VARIABLE);
	//status = H5Tset_cset(type, H5T_CSET_UTF8);
	
	dims[0] = nvarnames;
	dataspace_id = H5Screate_simple(1, dims, NULL);

	dataset_id = H5Dcreate(loc_id, VARIATE_NAMES, type, dataspace_id,
			       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	status = H5Dwrite(dataset_id, type, H5S_ALL, H5S_ALL, H5P_DEFAULT,
			  varnames);
	status = H5Dclose(dataset_id);
	status = H5Sclose(dataspace_id);

	free_varnames(varnames, nvarnames);

	return status;
}


static herr_t output_recv_fit(hid_t file_id, const struct recv_fit *fit, int has_resid)
{
	herr_t status = 0;
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_model *m = ll->model;
	struct frame *f = recv_model_frame(m);
	struct history *h = frame_history(f);

	hsize_t dim = (hsize_t)recv_model_dim(m);
	hsize_t nc = (hsize_t)constr_count(&constr);

	// intervals
	{
		hsize_t nintvl = (hsize_t)history_interval_count(h);
		const double *intvls = history_intervals(h);
		status = H5LTmake_dataset(file_id, INTERVALS, 1, &nintvl, H5T_NATIVE_DOUBLE, intvls);
	}

	// variate_names
	status = output_variate_names(file_id, f);

	// coefficients
	{
		const struct recv_coefs *coefs = recv_fit_coefs(fit);
		status = H5LTmake_dataset(file_id, COEFFICIENTS, 1, &dim, H5T_NATIVE_DOUBLE, coefs->all);
	}

	// constraints, constraint_values
	{
		hsize_t dims[2] = { dim, nc };
		const double *wts = constr_all_wts(&constr);
		const double *vals = constr_all_vals(&constr);

		status = H5LTmake_dataset(file_id, CONSTRAINTS, 2, dims, H5T_NATIVE_DOUBLE, wts);
		status = H5LTmake_dataset(file_id, CONSTRAINT_VALUES, 1, &nc, H5T_NATIVE_DOUBLE, vals);
	}

	// duals
	{
		const double *duals = recv_fit_duals(fit);
		status = H5LTmake_dataset(file_id, DUALS, 1, &nc, H5T_NATIVE_DOUBLE, duals);
	}

	// count
	{
		hsize_t count = recv_loglik_count(ll);
		status = H5LTmake_dataset(file_id, COUNT, 0, NULL, H5T_NATIVE_HSIZE, &count);
	}

	// score
	{
		struct recv_coefs score;
		recv_coefs_init(&score, f);
		memset(score.all, 0, dim * sizeof(double));
		recv_loglik_axpy_score(1.0, ll, &score);
		status = H5LTmake_dataset(file_id, SCORE, 1, &dim, H5T_NATIVE_DOUBLE, score.all);
		recv_coefs_deinit(&score);
	}

	// information
	{
		hsize_t dim2 = dim * (dim + 1) / 2;
		const char *uplo = (MLOGIT_COV_UPLO == BLAS_LOWER) ? "LOWER" : "UPPER";
		double *imat = xcalloc((size_t)dim2, sizeof(double));

		recv_loglik_axpy_imat(1.0, ll, imat);
		status = H5LTmake_dataset(file_id, INFORMATION, 1, &dim2, H5T_NATIVE_DOUBLE, imat);
		status = H5LTset_attribute_string(file_id, INFORMATION, "uplo", uplo);

		free(imat);
	}

	// rank, df_residual, df_null
	{
		hsize_t rank = dim - nc;
		hsize_t ntot = recv_loglik_count(ll);
		double dfresid = (double)(ntot - rank);
		double dfnull = (double)ntot;

		status = H5LTmake_dataset(file_id, RANK, 0, NULL, H5T_NATIVE_HSIZE, &rank);
		status = H5LTmake_dataset(file_id, DF_RESIDUAL, 0, NULL, H5T_NATIVE_DOUBLE, &dfresid);
		status = H5LTmake_dataset(file_id, DF_NULL, 0, NULL, H5T_NATIVE_DOUBLE, &dfnull);
	}

	// deviance, null_deviance
	{
		double dev = recv_fit_dev(fit);
		double dev0 = recv_fit_dev0(fit);
		
		status = H5LTmake_dataset(file_id, DEVIANCE, 0, NULL, H5T_NATIVE_DOUBLE, &dev);
		status = H5LTmake_dataset(file_id, NULL_DEVIANCE, 0, NULL, H5T_NATIVE_DOUBLE, &dev0);
	}


	// fitted_counts, observed_counts (residuals)
	if (has_resid) {
		const struct recv_coefs *coefs = recv_fit_coefs(fit);
		hsize_t dims[2];
		struct recv_resid resid;

		dims[0] = frame_send_count(f);
		dims[1] = frame_recv_count(f);

		history_clear(h);
		recv_resid_init(&resid, f, fit->ymsgs, coefs);

		status = H5LTmake_dataset(file_id, FITTED_COUNTS, 2, dims, H5T_NATIVE_DOUBLE, resid.fit.dyad);
		status = H5LTmake_dataset(file_id, OBSERVED_COUNTS, 2, dims, H5T_NATIVE_DOUBLE, resid.obs.dyad);
		recv_resid_deinit(&resid);
	}



	// effects
	// qr
	// linear.predictors (eta)
	// aic
	// iter
	// weights = wt
	// prior.weights = weights
	// converged = conv
	// boundary = boundary
	// call
	// formula
	// terms
	// data
	// offset
	// control
	// method
	// contrasts
	// xlevels

	return status;
}


static herr_t output(const char *output, const struct recv_fit *fit, int resid)
{
	hid_t file_id;
	herr_t status;

	// Create a new file using default properties.
	file_id = H5Fcreate(output, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	status = output_recv_fit(file_id, fit, resid);

	// Terminate access to the file.
	status = H5Fclose(file_id);
	
	return status;
}
*/

static int do_fit(struct recv_fit *fit, const struct recv_params *params0, const double *duals0)
{
	size_t maxit = 300;
	size_t report = 1;
	int trace = 1;
	enum recv_fit_task task;
	size_t it = 0;
	int err = 0;

	for (it = 0, task = recv_fit_start(fit, params0, duals0);
	     it < maxit && task == RECV_FIT_STEP;
	     it++, task = recv_fit_advance(fit)) {
		if (trace && it % report == 0 && task != RECV_FIT_CONV) {
			// const struct recv_loglik *ll = recv_fit_loglik(&fit);
			// const struct recv_loglik_info *info = recv_loglik_info(ll);
			// size_t n = info->nrecv;
			double dev = recv_fit_dev(fit);
			double nscore = recv_fit_score_norm(fit);
			double step = recv_fit_step_size(fit);

			fprintf(stderr, "iter %zu deviance = %.2f; |score| = %.16f; step = %.16f\n",
				it, dev, nscore, step);
		}
	}


	if (task != RECV_FIT_CONV) {
		fprintf(stderr, "ERROR: %s\n", "YOWZA!!"); //recv_fit_errmsg(fit));
		err = -1;
	}

	return err;
}


/*
static void init_params(struct recv_params *params,
			const struct frame *f,
			const struct constr *c,
			const char *filename)
{
	hid_t file_id;
	herr_t status;
	int rank;
	hsize_t dims;
	H5T_class_t type_class;
	size_t type_size;

	file_id = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0)
		exit(EXIT_FAILURE);

	recv_params_init(params, f, c);

	// coefficients
	{
		status = H5LTfind_dataset(file_id, COEFFICIENTS);
		if (!status)
			exit(EXIT_FAILURE);

		status = H5LTget_dataset_ndims(file_id, COEFFICIENTS, &rank);
		if (status < 0 || rank != 1)
			exit(EXIT_FAILURE);

		status = H5LTget_dataset_info(file_id, COEFFICIENTS, &dims, &type_class, &type_size);
		if (status < 0 || dims != params->coefs.dim)
			exit(EXIT_FAILURE);

		status = H5LTread_dataset(file_id, COEFFICIENTS, H5T_NATIVE_DOUBLE, params->coefs.all);
		if (status < 0)
			exit(EXIT_FAILURE);
	}

	// duals
	{
		status = H5LTfind_dataset(file_id, DUALS);
		if (!status)
			exit(EXIT_FAILURE);

		status = H5LTget_dataset_ndims(file_id, DUALS, &rank);
		if (status < 0 || rank != 1)
			exit(EXIT_FAILURE);

		status = H5LTget_dataset_info(file_id, DUALS, &dims, &type_class, &type_size);
		if (status < 0 || dims != constr_count(c))
			exit(EXIT_FAILURE);

		status = H5LTread_dataset(file_id, DUALS, H5T_NATIVE_DOUBLE, params->duals);
		if (status < 0)
			exit(EXIT_FAILURE);
	}

	status = H5Fclose(file_id);
	if (status < 0)
		exit(EXIT_FAILURE);
}

*/

/*
struct options {
	char *startfile;
	char *output;
	int boot;
	int resid;
	int32_t seed;
};

static struct options parse_options(int argc, char **argv)
{
	int c;
	static struct option long_options[] = {
		*//* These optoins set a flag. *//*
		{ "resid", no_argument, 0, 'r' },

		*//* These options don't set a flag *//*
		{ "output", required_argument, 0, 'o' },
		{ "start", required_argument, 0, 's' },
		{ "boot", optional_argument, 0, 'b' },
		{ NULL, 0, NULL, 0 }
	};
	int option_index = 0;

	struct options opts;
	opts.output = "output.h5";
	opts.startfile = NULL;
	opts.boot = 0;
	opts.resid = 0;
	opts.seed = 0;

	while ((c = getopt_long(argc, argv, "ro:s:b::", long_options, &option_index)) != -1) {
		switch (c) {
		case 's':
			opts.startfile = optarg;
			fprintf(stderr, "start file: '%s'\n", opts.startfile);
			fflush(stderr);
			break;
		case 'b':
			opts.boot = 1;
			if (optarg) {
				opts.seed = (int32_t)strtol(optarg, NULL, 10);
			}
			fprintf(stderr, "bootstrap seed: '%zu'\n", (size_t)opts.seed);
			fflush(stderr);
			break;
		case 'r':
			opts.resid = 1;
			break;
		case 'o':
			opts.output = optarg;
			break;
		default:
			exit(1);
		}
	}

	if (opts.boot && !opts.startfile) {
		fprintf(stderr, "Must specify start file for bootstrap\n");
		exit(1);
	}

	return opts;
}
*/

int main(int argc, char **argv)
{
	/* struct options opts = parse_options(argc, argv); */
	struct history h;
	struct design r;
	struct design2 d;
	struct constr c;
	int exclude_loops = 1;
	struct recv_fit fit;
	const struct message *msgs;
	size_t nmsg;
	size_t ncextra;
	int err;

	err = setup_history(&h);
	if (err)
		goto fail_history;

	err = setup_designs(&r, &d, &h);
	if (err)
		goto fail_designs;

	setup_constr(&c, &r, &d);

	/*
	// setup initial parameters
	struct recv_params params0;
	int has_params0 = 0;
	if (opts.startfile) {
		has_params0 = 1;
		init_params(&params0, &frame, &constr, opts.startfile);
	}
	struct recv_params *pparams0 = has_params0 ? &params0 : NULL;
	struct recv_coefs *pcoefs0 = has_params0 ? &params0.coefs : NULL;

	// setup bootstrap (if required)
	struct recv_boot boot;
	if (opts.boot) {
		dsfmt_t dsfmt;
		dsfmt_init_gen_rand(&dsfmt, opts.seed);

		recv_boot_init(&boot, &frame, &enron_messages, pcoefs0, &dsfmt);
		ymsgs = &boot.messages;
	}
	*/

	history_get_messages(&h, &msgs, &nmsg);
	recv_fit_init(&fit, &r, &d, exclude_loops, msgs, nmsg, &c, NULL);

	ncextra = recv_fit_extra_constr_count(&fit);
	if (ncextra)
		fprintf(stderr, "Adding %zd %s to make parameters identifiable\n",
			ncextra, ncextra == 1 ? "constraint" : "constraints");

	// fit the model
	err = do_fit(&fit, NULL, NULL);

	// output the results
	//output(opts.output, &fit, opts.resid);

	//if (opts.boot)
	//	recv_boot_deinit(&boot);
	//if (has_params0)
	//	recv_params_deinit(&params0);

	recv_fit_deinit(&fit);
	constr_deinit(&c);
	design2_deinit(&d);
	design_deinit(&r);
fail_designs:
	history_deinit(&h);
fail_history:
	return err;
}
