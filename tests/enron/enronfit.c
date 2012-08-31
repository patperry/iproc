#include "port.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
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
	size_t terms = 2; // second-order interactions only
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
		// 58982400.00
	};
	size_t nintvls = sizeof(intvls) / sizeof(intvls[0]);
	//int has_effects = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, nintvls);

	/* send design */
	struct design *s = frame_send_design(&frame);
	design_add_traits(s, trait_names, traits, ntrait);
       
	/* recv design */
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);
	
	/* dyad design */
	struct design2 *d = frame_dyad_design(&frame);

	design2_add_tvar(d, "IRecv", DYAD_VAR_IRECV);
	design2_add_tvar(d, "NRecv", DYAD_VAR_NRECV);
	design2_add_tvar(d, "ISend", DYAD_VAR_ISEND);
	design2_add_tvar(d, "NSend", DYAD_VAR_NSEND);
	//design2_add_tvar(d, "IRecv2", DYAD_VAR_IRECV2);
	//design2_add_tvar(d, "NRecv2", DYAD_VAR_NRECV2);
	//design2_add_tvar(d, "ISend2", DYAD_VAR_ISEND2);
	//design2_add_tvar(d, "NSend2", DYAD_VAR_NSEND2);
	//design2_add_tvar(d, "ISib", DYAD_VAR_ISIB);
	//design2_add_tvar(d, "NSib", DYAD_VAR_NSIB);
	//design2_add_tvar(d, "ICosib", DYAD_VAR_ICOSIB);
	//design2_add_tvar(d, "NCosib", DYAD_VAR_NCOSIB);


	char buf[1024];
	size_t i, j;
	for (i = 0; i < ntrait; i++) {
		for (j = 0; j < ntrait; j++) {
			sprintf(buf, "%s*%s", trait_names[i], trait_names[j]);
			design2_add_kron(d, buf, design_var(s, trait_names[i]), design_var(r, trait_names[j]));
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


#define YG(gen) \
	do { \
		if ((err = gen) != yajl_gen_status_ok) \
			return err; \
	} while (0)

#define YSTR(str) ((const unsigned char *)(str)), (strlen(str))

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
#define OBSERVED_COUNTS		"observed_counts"
#define EXPECTED_COUNTS		"expected_counts"


yajl_gen_status yajl_gen_recv_fit(yajl_gen hand, const struct recv_fit *fit,
				  int has_resid)
{
	assert(fit);
	yajl_gen_status err = yajl_gen_status_ok;
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_model *m = ll->model;
	struct frame *f = recv_model_frame(m);

	size_t ne = recv_fit_constr_count(fit);
	size_t dim = recv_model_dim(m);
	size_t i, n;

	YG(yajl_gen_map_open(hand));
	{
		// intervals
		YG(yajl_gen_string(hand, YSTR(INTERVALS)));
		const double *intvls = frame_intervals(f);
		n = frame_interval_count(f);
		YG(yajl_gen_array_open(hand));
		for (i = 0; i < n; i++) {
			YG(yajl_gen_ieee754(hand, intvls[i]));
		}
		YG(yajl_gen_array_close(hand));

		// variate_names
		//YG(yajl_gen_string(hand, YSTR(VARIATE_NAMES)));
		//YG(yajl_gen_array_open(hand));
		//n = design_dim(d);
		//for (i = 0; i < n; i++) {
		//		YG(yajl_gen_string(hand, YSTR(design_name(d, i))));
		//}
		//YG(yajl_gen_array_close(hand));

		// coefficients
		YG(yajl_gen_string(hand, YSTR(COEFFICIENTS)));
		const struct recv_coefs *coefs = recv_fit_coefs(fit);
		YG(yajl_gen_vector(hand, dim, coefs->all));

		// constraint names
		YG(yajl_gen_string(hand, YSTR(CONSTRAINT_NAMES)));
		YG(yajl_gen_array_open(hand));
		{
			n = recv_fit_constr_count(fit);

			for (i = 0; i < n; i++) {
				const char *name;
				recv_fit_get_constr(fit, i, NULL, NULL, NULL, NULL, &name);
				YG(yajl_gen_string(hand, YSTR(name)));
			}
		}
		YG(yajl_gen_array_close(hand));

		// constraints
		YG(yajl_gen_string(hand, YSTR(CONSTRAINTS)));
		YG(yajl_gen_array_open(hand));
		{
			n = recv_fit_constr_count(fit);

			for (i = 0; i < n; i++) {
				const double *wts;
				const size_t *ind;
				size_t iz, nz = 0;
				recv_fit_get_constr(fit, i, &wts, &ind, &nz, NULL, NULL);

				YG(yajl_gen_map_open(hand));
				YG(yajl_gen_string(hand, YSTR("dim")));
				YG(yajl_gen_integer(hand, dim));

				YG(yajl_gen_string(hand, YSTR("count")));
				YG(yajl_gen_integer(hand, nz));

				YG(yajl_gen_string(hand, YSTR("pattern")));
				YG(yajl_gen_array_open(hand));
				for (iz = 0; iz < nz; iz++) {
					YG(yajl_gen_integer(hand, ind[iz] + 1));
				}
				YG(yajl_gen_array_close(hand));

				YG(yajl_gen_string(hand, YSTR("data")));
				YG(yajl_gen_array_open(hand));
				for (iz = 0; iz < nz; iz++) {
					YG(yajl_gen_ieee754(hand, wts[iz]));
				}
				YG(yajl_gen_array_close(hand));
				YG(yajl_gen_map_close(hand));
			}
		}
		YG(yajl_gen_array_close(hand));

		// constraint values
		YG(yajl_gen_string(hand, YSTR(CONSTRAINT_VALUES)));
		YG(yajl_gen_array_open(hand));
		{
			n = recv_fit_constr_count(fit);

			for (i = 0; i < n; i++) {
				double value;
				recv_fit_get_constr(fit, i, NULL, NULL, NULL, &value, NULL);
				YG(yajl_gen_ieee754(hand, value));
			}
		}
		YG(yajl_gen_array_close(hand));

		// duals
		YG(yajl_gen_string(hand, YSTR(DUALS)));
		const double *duals = recv_fit_duals(fit);
		YG(yajl_gen_vector(hand, ne, duals));

		// count
		YG(yajl_gen_string(hand, YSTR(COUNT)));
		size_t count = recv_loglik_count(ll);
		YG(yajl_gen_integer(hand, count));

		// score
		struct recv_coefs score;
		
		recv_coefs_init(&score, f);
		memset(score.all, 0, dim * sizeof(double));
		recv_loglik_axpy_score(1.0, ll, &score);
		YG(yajl_gen_string(hand, YSTR(SCORE)));
		YG(yajl_gen_vector(hand, dim, score.all));

		recv_coefs_deinit(&score);

		// information
		double *imat = xcalloc(dim * (dim + 1) / 2, sizeof(double));
		recv_loglik_axpy_imat(1.0, ll, imat);

		YG(yajl_gen_string(hand, YSTR(INFORMATION)));
		YG(yajl_gen_vector(hand, dim * (dim + 1) / 2, imat));

		free(imat);

		// rank
		YG(yajl_gen_string(hand, YSTR(RANK)));
		size_t rank = dim - ne;
		YG(yajl_gen_integer(hand, rank));

		// deviance
		YG(yajl_gen_string(hand, YSTR(DEVIANCE)));
		double dev = recv_fit_dev(fit);
		YG(yajl_gen_ieee754(hand, dev));

		// df.residual
		YG(yajl_gen_string(hand, YSTR(DF_RESIDUAL)));
		size_t ntot = recv_loglik_count(ll);
		YG(yajl_gen_integer(hand, ntot - rank));

		// null.deviance
		YG(yajl_gen_string(hand, YSTR(NULL_DEVIANCE)));
		YG(yajl_gen_ieee754(hand, recv_fit_dev0(fit)));

		// df.null
		YG(yajl_gen_string(hand, YSTR(DF_NULL)));
		YG(yajl_gen_integer(hand, ntot));

		// residuals
		// fitted.values
		// y = y
		if (has_resid) {
			struct recv_resid resid;

			frame_clear(f);
			recv_resid_init(&resid, f, fit->ymsgs, coefs);
			YG(yajl_gen_string(hand, YSTR(OBSERVED_COUNTS)));
			YG(yajl_gen_matrix(hand, nsend, nrecv, resid.obs.dyad));

			YG(yajl_gen_string(hand, YSTR(EXPECTED_COUNTS)));
			YG(yajl_gen_matrix(hand, nsend, nrecv, resid.exp.dyad));

			recv_resid_deinit(&resid);
		}
	}
	YG(yajl_gen_map_close(hand));

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


	return err;
}


static void output(const struct recv_fit *fit, int resid)
{
	/* generator config */
	yajl_gen g = yajl_gen_alloc(NULL);
	yajl_gen_config(g, yajl_gen_beautify, 1);
	yajl_gen_recv_fit(g, fit, resid);


	/* output */
	const unsigned char * buf;
	size_t len;
	yajl_gen_get_buf(g, &buf, &len);
	fwrite(buf, 1, len, stdout);
	yajl_gen_clear(g);
	yajl_gen_free(g);
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
	int err = 1;
	FILE *fp = fopen(filename, "rb");
	if (!fp)
		goto out;

	fseek(fp, 0L, SEEK_END);
	size_t sz = ftell(fp);
	fseek(fp, 0L, SEEK_SET);

	char *filebuf = xmalloc(sz + 1);
	size_t rd = fread(filebuf, 1, sz, fp);

	// file read error handling
	if (rd == 0 && !feof(fp)) {
		fprintf(stderr, "error encountered on file read\n");
		goto cleanup_file;
	}
	assert(rd == sz);
	filebuf[sz] = '\0';

	char errbuf[1024];
	errbuf[0] = '\0';

	yajl_val node = yajl_tree_parse(filebuf, errbuf, sizeof(errbuf));

	// parse error handling
	if (node == NULL) {
		fprintf(stderr, "parse_error: ");
		if (strlen(errbuf)) fprintf(stderr, " %s", errbuf);
		else fprintf(stderr, "unknown error");
		fprintf(stderr, "\n");
		goto cleanup_yajl;
	}

	// ... and extract a nested value from the config file
	const char * coefs_path[] = { COEFFICIENTS, NULL };
	yajl_val vcoefs = yajl_tree_get(node, coefs_path,  yajl_t_array);

	if (!vcoefs
	    || YAJL_GET_ARRAY(vcoefs)->len != recv_fit_coefs_count(fit))
		goto cleanup_yajl;

	size_t i, dim = YAJL_GET_ARRAY(vcoefs)->len;
	double *coefs = xcalloc(dim, sizeof(*coefs));
	yajl_val *item = YAJL_GET_ARRAY(vcoefs)->values;

	for (i = 0; i < dim; i++) {
		if (!YAJL_IS_NUMBER(item[i])) {
			err = 1;
			goto cleanup_coefs;
		}
		coefs[i] = YAJL_GET_DOUBLE(item[i]);
	}
	
	const char * duals_path[] = { DUALS, NULL };
	yajl_val vduals = yajl_tree_get(node, duals_path,  yajl_t_array);

	if (!vduals
	    || YAJL_GET_ARRAY(vduals)->len != recv_fit_constr_count(fit))
		goto cleanup_coefs;

	size_t nduals = YAJL_GET_ARRAY(vduals)->len;
	double *duals = xcalloc(nduals, sizeof(*duals));
	item = YAJL_GET_ARRAY(vduals)->values;

	for (i = 0; i < nduals; i++) {
		if (!YAJL_IS_NUMBER(item[i])) {
			err = 1;
			goto cleanup_duals;
		}
		duals[i] = YAJL_GET_DOUBLE(item[i]);
	}

	recv_fit_params_init(params, fit);
	memcpy(params->coefs.all, coefs, dim * sizeof(double));
	memcpy(params->duals, duals, nduals * sizeof(double));
	err = 0;

cleanup_duals:
	free(duals);
cleanup_coefs:
	free(coefs);
cleanup_yajl:
	yajl_tree_free(node);
	free(filebuf);
cleanup_file:
	fclose(fp);
out:
	if (err)
		exit(err);
}


struct options {
	char *startfile;
	int boot;
	int resid;
	int32_t seed;
};

static struct options parse_options(int argc, char **argv)
{
	int c;
	static struct option long_options[] = {
		{ "start", required_argument, 0, 's' },
		{ "boot", optional_argument, 0, 'b' },
		{ "resid", no_argument, 0, 'r' },
		{ NULL, 0, NULL, 0 }
	};
	int option_index = 0;

	struct options opts;
	opts.startfile = NULL;
	opts.boot = 0;
	opts.resid = 0;
	opts.seed = 0;


	while ((c = getopt_long(argc, argv, "s:b:", long_options, &option_index)) != -1) {
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
	output(&fit, opts.resid);

	if (opts.boot)
		recv_boot_deinit(&boot);
	if (has_params0)
		recv_fit_params_deinit(&params0);
	recv_fit_deinit(&fit);
	messages_deinit(&enron_messages);
	teardown_frame();

	return err;
}
