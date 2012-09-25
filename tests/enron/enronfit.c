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
#include "yajl/yajl_tree.h"
#define DSFMT_MEXP 19937
#define DSFMT_DO_NOT_USE_OLD_NAMES
#include "dSFMT/dSFMT.h"

#include "ieee754.h"
#include "json.h"
#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"
#include "recv_boot.h"
#include "recv_model.h"
#include "recv_loglik.h"
#include "recv_fit.h"
#include "recv_resid.h"


static size_t nsend;
static size_t nrecv;
static size_t ntrait;
static double *traits;
static const char * const *trait_names;
static int has_loops;

static struct frame frame;


static void setup_frame(void) {
	char namebuf[1024];
	size_t i, j;

	size_t terms = 1; // first-order interactions only
	enron_employees_init(&nsend, &traits, &ntrait, &trait_names, terms);
	nrecv = nsend;
	has_loops = 0;

	double intvls[] = {
		// 450.00,
		//900.00,    // 15 min
		1800.00,   // 30 min
			   //3600.00,   //  1 hr
		7200.00,   //  2 hr
			   //14400.00,  //  4 hr
		28800.00,  //  8 hr
			   //57600.00,  // 16 hr
		115200.00, // 32 hr
			   // 230400.00, // 2.66 day 
		460800.00, // 5.33 day
			   // 921600.00, // 10.66 day
		1843200.00, // 21.33 day
			    // 3686400.00, // 42.66 day
			    // 7372800.00, // 85.33 day
			    // 14745600.00, // 170.66 day
			    // 29491200.00 // 341.33 day
			    // 58982400.00 // 682.66 day
		INFINITY
	};
	size_t nintvls = sizeof(intvls) / sizeof(intvls[0]);
	//int has_effects = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, nintvls);

	/* send design */
	struct design *s = frame_send_design(&frame);
	design_add_traits(s, trait_names, traits, ntrait);
	design_add_prod(s, "Leg:Jun", design_var(s, "Leg"), design_var(s, "Jun"));
	design_add_prod(s, "Trad:Jun", design_var(s, "Trad"), design_var(s, "Jun"));
	design_add_prod(s, "Leg:Fem", design_var(s, "Leg"), design_var(s, "Fem"));
	design_add_prod(s, "Trad:Fem", design_var(s, "Trad"), design_var(s, "Fem"));
	design_add_prod(s, "Jun:Fem", design_var(s, "Jun"), design_var(s, "Fem"));
       
	/* recv design */
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);
	design_add_prod(r, "Leg:Jun", design_var(r, "Leg"), design_var(r, "Jun"));
	design_add_prod(r, "Trad:Jun", design_var(r, "Trad"), design_var(r, "Jun"));
	design_add_prod(r, "Leg:Fem", design_var(r, "Leg"), design_var(r, "Fem"));
	design_add_prod(r, "Trad:Fem", design_var(r, "Trad"), design_var(r, "Fem"));
	design_add_prod(r, "Jun:Fem", design_var(r, "Jun"), design_var(r, "Fem"));

	design_add_tvar(r, "ISendTot", VAR_ISENDTOT);
	design_add_tvar(r, "IRecvTot", VAR_IRECVTOT);
	design_add_prod(r, "ISendTot:IRecvTot", design_var(r, "ISendTot"), design_var(r, "IRecvTot"));

	design_add_tvar(r, "NRecvTot", VAR_NRECVTOT);
	design_add_tvar(r, "NSendTot", VAR_NSENDTOT);

	/* dyad design */
	struct design2 *d = frame_dyad_design(&frame);

	design2_add_tvar(d, "ISend", VAR2_ISEND);
	design2_add_tvar(d, "IRecv", VAR2_IRECV);
	design2_add_prod(d, "ISend:IRecv", design2_var(d, "ISend"), design2_var(d, "IRecv"));

	design2_add_tvar(d, "NSend", VAR2_NSEND);
	design2_add_tvar(d, "NRecv", VAR2_NRECV);

	//design2_add_tvar(d, "ISend2", VAR2_ISEND2);
	//design2_add_tvar(d, "IRecv2", VAR2_IRECV2);
	//design2_add_tvar(d, "ISib", VAR2_ISIB);
	//design2_add_tvar(d, "ICosib", VAR2_ICOSIB);

	//design2_add_tvar(d, "NSend2", VAR2_NSEND2);
	//design2_add_tvar(d, "NRecv2", VAR2_NRECV2);
	//design2_add_tvar(d, "NSib", VAR2_NSIB);
	//design2_add_tvar(d, "NCosib", VAR2_NCOSIB);


	for (i = 0; i < design_trait_count(s); i++) {
		const struct var *vi = design_trait_var(s, i);
		for (j = 0; j < design_trait_count(r); j++) {
			const struct var *vj = design_trait_var(r, j);
			sprintf(namebuf, "%s*%s", var_name(vi), var_name(vj));
			design2_add_kron(d, namebuf, vi, vj);
		}
	}



	// recv_model_init(&m, &frame);
	//recv_model_add_inter(&m, design_var(r, "Female"), design_var(d, "NRecv"));
	
	/* dyad design */
	//struct dyad_design *d = frame_dyad_design(&frame);
	//dyad_design_add_interact(d, design_ix(s, "Legal"), design_ix(r, "Legal"));
}

static void add_constraints(struct recv_fit *fit)
{
	(void)fit;
	/* add constraints to make the model identifiable (TODO/not implemented) */
	//size_t nadd = recv_fit_add_constr_identify(fit);
	//if (nadd > 0)
	//	fprintf(stderr, "Adding %d constraints to make parameters identifiable\n", nadd);
}

static void teardown_frame(void)
{
	frame_deinit(&frame);
	free(traits);
}


#define INTERVALS		"intervals"
#define VARIATE_NAMES		"variate_names"
#define COHORT_NAMES		"cohort_names"
#define COEFFICIENTS		"coefficients"
#define COUNT			"count"
#define SCORE			"score"
#define INFORMATION		"information"
#define RANK			"rank"
#define CONSTRAINT_NAMES	"constraint_names"
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
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_dim(r);
	size_t dimd = design2_dim(d);
	size_t i, dim = dimr + dimd;

	const char **res = xmalloc(dim * sizeof(*res));
	for (i = 0; i < dim; i++) {
		res[i] = "";
	}

	*pn = dim;

	return res;
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

	free(varnames);

	return status;
}


static herr_t output_recv_fit(hid_t file_id, const struct recv_fit *fit, int has_resid)
{
	herr_t status = 0;
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_model *m = ll->model;
	struct frame *f = recv_model_frame(m);

	hsize_t ne = (hsize_t)recv_fit_constr_count(fit);
	hsize_t dim = (hsize_t)recv_model_dim(m);
	hsize_t one = 1;
	//size_t i, n;

	// intervals
	{
		hsize_t nintvl = (hsize_t)frame_interval_count(f);
		const double *intvls = frame_intervals(f);
		status = H5LTmake_dataset(file_id, INTERVALS, 1, &nintvl, H5T_NATIVE_DOUBLE, intvls);
	}

	// variate_names
	status = output_variate_names(file_id, f);

	// coefficients
	{
		const struct recv_coefs *coefs = recv_fit_coefs(fit);
		status = H5LTmake_dataset(file_id, COEFFICIENTS, 1, &dim, H5T_NATIVE_DOUBLE, coefs->all);
	}

	// constraint_names
	// constraints
	// constraint_values

	// duals
	{
		const double *duals = recv_fit_duals(fit);
		status = H5LTmake_dataset(file_id, DUALS, 1, &ne, H5T_NATIVE_DOUBLE, duals);
	}

	// count
	{
		hsize_t count = recv_loglik_count(ll);
		status = H5LTmake_dataset(file_id, COUNT, 1, &one, H5T_NATIVE_HSIZE, &count);
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
		hsize_t rank = dim - ne;
		hsize_t ntot = recv_loglik_count(ll);
		double dfresid = (double)(ntot - rank);
		double dfnull = (double)ntot;

		status = H5LTmake_dataset(file_id, RANK, 1, &one, H5T_NATIVE_HSIZE, &rank);
		status = H5LTmake_dataset(file_id, DF_RESIDUAL, 1, &one, H5T_NATIVE_DOUBLE, &dfresid);
		status = H5LTmake_dataset(file_id, DF_NULL, 1, &one, H5T_NATIVE_DOUBLE, &dfnull);
	}

	// deviance, null_deviance
	{
		double dev = recv_fit_dev(fit);
		double dev0 = recv_fit_dev0(fit);
		
		status = H5LTmake_dataset(file_id, DEVIANCE, 1, &one, H5T_NATIVE_DOUBLE, &dev);
		status = H5LTmake_dataset(file_id, NULL_DEVIANCE, 1, &one, H5T_NATIVE_DOUBLE, &dev0);
	}


	// fitted_counts, observed_counts (residuals)
	if (has_resid) {
		const struct recv_coefs *coefs = recv_fit_coefs(fit);
		hsize_t dims[2];
		struct recv_resid resid;

		dims[0] = frame_send_count(f);
		dims[1] = frame_recv_count(f);


		frame_clear(f);
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

	/* Create a new file using default properties. */
	file_id = H5Fcreate(output, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	status = output_recv_fit(file_id, fit, resid);

	/* Terminate access to the file. */
	status = H5Fclose(file_id);
	
	return status;
}

static int do_fit(struct recv_fit *fit, const struct recv_fit_params *params0)
{
	size_t maxit = 30;
	size_t report = 1;
	int trace = 1;


	enum recv_fit_task task;
	size_t it = 0;

	//struct matrix coefs0;
	//matrix_init(&coefs0, dim, nc);
	//matrix_fill(&coefs0, 0.0001);

	for (it = 0, task = recv_fit_start(fit, params0);
	     it < maxit && task == RECV_FIT_STEP;
	     it++, task = recv_fit_advance(fit)) {
		if (trace && it % report == 0 && task != RECV_FIT_CONV) {
			// const struct recv_loglik *ll = recv_fit_loglik(&fit);
			// const struct recv_loglik_info *info = recv_loglik_info(ll);
			// size_t n = info->nrecv;
			// double dev = n * info->dev;
			// const struct vector *score = &info->score;
			// double ngrad = vector_max_abs(score);
			double step = recv_fit_step(fit);
			double ngrad = recv_fit_grad_norm2(fit);
			fprintf(stderr, "iter %zu; |grad| = %.16f; step = %.16f\n",
				it, ngrad, step);

			// fprintf(stderr, "iter %zu deviance = %.2f; |grad| = %.16f; step = %.16f\n",
			//	it, dev, ngrad, step);
		}
	}

	//matrix_deinit(&coefs0);

	int err = 0;

	if (task != RECV_FIT_CONV) {
		fprintf(stderr, "ERROR: %s\n", recv_fit_errmsg(fit));
		err = -1;
	}

	return err;
}


static void init_params(struct recv_fit_params *params,
			const struct recv_fit *fit,
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

	recv_fit_params_init(params, fit);

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
		if (status < 0 || dims != recv_fit_constr_count(fit))
			exit(EXIT_FAILURE);

		status = H5LTread_dataset(file_id, DUALS, H5T_NATIVE_DOUBLE, params->duals);
		if (status < 0)
			exit(EXIT_FAILURE);
	}

	status = H5Fclose(file_id);
	if (status < 0)
		exit(EXIT_FAILURE);
}


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
		/* These optoins set a flag. */
		{ "resid", no_argument, 0, 'r' },

		/* These options don't set a flag */
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


int main(int argc, char **argv)
{
	int err = 0;

	struct options opts = parse_options(argc, argv);

	// setup frame
	setup_frame();

	// setup messages
	struct messages enron_messages;
	struct messages *xmsgs = &enron_messages;
	struct messages *ymsgs = &enron_messages;
	enron_messages_init(&enron_messages, 5);

	// setup fit
	struct recv_fit fit;
	struct recv_fit_ctrl ctrl = RECV_FIT_CTRL0;
	recv_fit_init(&fit, &frame, xmsgs, ymsgs, &ctrl);
	add_constraints(&fit);

	// setup initial parameters
	struct recv_fit_params params0;
	int has_params0 = 0;
	if (opts.startfile) {
		has_params0 = 1;
		init_params(&params0, &fit, opts.startfile);
	}
	struct recv_fit_params *pparams0 = has_params0 ? &params0 : NULL;
	struct recv_coefs *pcoefs0 = has_params0 ? &params0.coefs : NULL;

	// setup bootstrap (if required)
	struct recv_boot boot;
	if (opts.boot) {
		dsfmt_t dsfmt;
		dsfmt_init_gen_rand(&dsfmt, opts.seed);

		recv_boot_init(&boot, &frame, &enron_messages, pcoefs0, &dsfmt);
		ymsgs = &boot.messages;
	}

	// fit the model
	err = do_fit(&fit, pparams0);

	// output the results
	output(opts.output, &fit, opts.resid);

	if (opts.boot)
		recv_boot_deinit(&boot);
	if (has_params0)
		recv_fit_params_deinit(&params0);
	recv_fit_deinit(&fit);
	messages_deinit(&enron_messages);
	teardown_frame();

	return err;
}
