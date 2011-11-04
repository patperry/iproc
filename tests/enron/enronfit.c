#include "port.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <getopt.h>
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
static size_t ncohort;
static size_t ntrait;
static size_t *cohorts;
static double *traits;
static const char * const *cohort_names;
static const char * const *trait_names;
static int has_loops;

static struct frame frame;


static void setup(void) {
	enron_employees_init(&nsend, &cohorts, &ncohort, &cohort_names,
			     &traits, &ntrait, &trait_names);
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
	int has_effects = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, nintvls);
	struct design *d = frame_recv_design(&frame);
        design_set_has_effects(d, has_effects);
	design_set_traits(d, traits, ntrait, trait_names);
	design_add_dvar(d, RECV_VAR_IRECV, NULL);
	design_add_dvar(d, RECV_VAR_NRECV, NULL);
	design_add_dvar(d, RECV_VAR_ISEND, NULL);
	design_add_dvar(d, RECV_VAR_NSEND, NULL);
	design_add_dvar(d, RECV_VAR_IRECV2, NULL);
	design_add_dvar(d, RECV_VAR_NRECV2, NULL);
	design_add_dvar(d, RECV_VAR_ISEND2, NULL);
	design_add_dvar(d, RECV_VAR_NSEND2, NULL);
	design_add_dvar(d, RECV_VAR_ISIB, NULL);
	design_add_dvar(d, RECV_VAR_NSIB, NULL);
	design_add_dvar(d, RECV_VAR_ICOSIB, NULL);
	design_add_dvar(d, RECV_VAR_NCOSIB, NULL);
}

static void add_constraints(struct recv_fit *fit)
{
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_model *m = ll->model;
	int dim = (int)recv_model_dim(m);
	size_t nc = recv_model_cohort_count(m);

	double ce[dim * nc];
	size_t ind[dim * nc];
	size_t nz;
	int i;
	char name[1024];


	for (i = 0; i < dim; i++) {

		// intercept
		if (!((0 <= i && i < 9) // 2nd-order interactions
		      || (11 <= i) // dynamic effects
		      )) {

			ind[0] = i + dim * ENRON_COHORT_OSM; ce[0] = 1.0;
			nz = 1;
			sprintf(name, "Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);
		}

		if (!((0 <= i && i < 9) // 2nd-order interactions
		      )) {
			ind[0] = i + dim * ENRON_COHORT_LSM; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_OSM; ce[1] = -1.0;
			nz = 2;
			sprintf(name, "L*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_TSM; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_OSM; ce[1] = -1.0;
			nz = 2;
			sprintf(name, "T*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_OJM; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_OSM; ce[1] = -1.0;
			nz = 2;
			sprintf(name, "J*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_OSF; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_OSM; ce[1] = -1.0;
			nz = 2;
			sprintf(name, "F*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			// 2-way interactions
			ind[0] = i + dim * ENRON_COHORT_LJM; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_LSM; ce[1] = -1.0;
			ind[2] = i + dim * ENRON_COHORT_OJM; ce[2] = -1.0;
			ind[3] = i + dim * ENRON_COHORT_OSM; ce[3] = +1.0;
			nz = 4;
			sprintf(name, "LJ*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_TJM; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_TSM; ce[1] = -1.0;
			ind[2] = i + dim * ENRON_COHORT_OJM; ce[2] = -1.0;
			ind[3] = i + dim * ENRON_COHORT_OSM; ce[3] = +1.0;
			nz = 4;
			sprintf(name, "TJ*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_LSF; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_LSM; ce[1] = -1.0;
			ind[2] = i + dim * ENRON_COHORT_OSF; ce[2] = -1.0;
			ind[3] = i + dim * ENRON_COHORT_OSM; ce[3] = +1.0;
			nz = 4;
			sprintf(name, "LF*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_TSF; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_TSM; ce[1] = -1.0;
			ind[2] = i + dim * ENRON_COHORT_OSF; ce[2] = -1.0;
			ind[3] = i + dim * ENRON_COHORT_OSM; ce[3] = +1.0;
			nz = 4;
			sprintf(name, "TF*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

			ind[0] = i + dim * ENRON_COHORT_OJF; ce[0] = +1.0;
			ind[1] = i + dim * ENRON_COHORT_OJM; ce[1] = -1.0;
			ind[2] = i + dim * ENRON_COHORT_OSF; ce[2] = -1.0;
			ind[3] = i + dim * ENRON_COHORT_OSM; ce[3] = +1.0;
			nz = 4;
			sprintf(name, "JF*Var%d", i + 1);
			recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);
		}

		// 3-way interactions
		ind[0] = i + dim * ENRON_COHORT_LJF; ce[0] = +1.0;
		ind[1] = i + dim * ENRON_COHORT_LJM; ce[1] = -1.0;
		ind[2] = i + dim * ENRON_COHORT_LSF; ce[2] = -1.0;
		ind[3] = i + dim * ENRON_COHORT_LSM; ce[3] = +1.0;
		ind[4] = i + dim * ENRON_COHORT_OJF; ce[4] = -1.0;
		ind[5] = i + dim * ENRON_COHORT_OJM; ce[5] = +1.0;
		ind[6] = i + dim * ENRON_COHORT_OSF; ce[6] = +1.0;
		ind[7] = i + dim * ENRON_COHORT_OSM; ce[7] = -1.0;
		nz = 8;
		sprintf(name, "LJF*Var%d", i + 1);
		recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);

		ind[0] = i + dim * ENRON_COHORT_TJF; ce[0] = +1.0;
		ind[1] = i + dim * ENRON_COHORT_TJM; ce[1] = -1.0;
		ind[2] = i + dim * ENRON_COHORT_TSF; ce[2] = -1.0;
		ind[3] = i + dim * ENRON_COHORT_TSM; ce[3] = +1.0;
		ind[4] = i + dim * ENRON_COHORT_OJF; ce[4] = -1.0;
		ind[5] = i + dim * ENRON_COHORT_OJM; ce[5] = +1.0;
		ind[6] = i + dim * ENRON_COHORT_OSF; ce[6] = +1.0;
		ind[7] = i + dim * ENRON_COHORT_OSM; ce[7] = -1.0;
		nz = 8;
		sprintf(name, "TJF*Var%d", i + 1);
		recv_fit_add_constr(fit, ce, ind, nz, 0.0, name);
	}

	/* add constraints to make the model identifiable (TODO/not implemented) */
	//size_t nadd = recv_fit_add_constr_identify(fit);
	//if (nadd > 0)
	//	fprintf(stderr, "Adding %d constraints to make parameters identifiable\n", nadd);
}

static void teardown(void)
{
	frame_deinit(&frame);
	free(traits);
	free(cohorts);
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

yajl_gen_status yajl_gen_recv_fit(yajl_gen hand, const struct recv_fit *fit)
{
	assert(fit);
	yajl_gen_status err = yajl_gen_status_ok;
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_model *m = ll->model;
	struct frame *f = recv_model_frame(m);
	const struct design *d = recv_model_design(m);

	size_t ne = recv_fit_constr_count(fit);
	size_t dim = recv_model_dim(m);
	size_t ic, nc = recv_model_cohort_count(m);
	size_t i, n;

	YG(yajl_gen_map_open(hand));
	{
		/* intervals */
		YG(yajl_gen_string(hand, YSTR(INTERVALS)));
		const double *intvls = frame_intervals(f);
		n = frame_interval_count(f);
		YG(yajl_gen_array_open(hand));
		for (i = 0; i < n; i++) {
			YG(yajl_gen_ieee754(hand, intvls[i]));
		}
		YG(yajl_gen_array_close(hand));

		/* variate_names */
		YG(yajl_gen_string(hand, YSTR(VARIATE_NAMES)));
		YG(yajl_gen_array_open(hand));
		const char * const *trait_names = design_trait_names(d);
		n = design_traits_dim(d);
		assert(design_traits_index(d) == 0);
		for (i = 0; i < n; i++) {
			YG(yajl_gen_string(hand, YSTR(trait_names[i])));
		}

		assert(design_dvars_index(d) == design_traits_dim(d));
		const struct design_var *rvs;
		size_t irv, nrv;
		design_get_dvars(d, &rvs, &nrv);
		for (irv = 0; irv < nrv; irv++) {
			n = rvs[irv].dim;
			for (i = 0; i < n; i++) {
				const char *name = rvs[irv].names[i];
				YG(yajl_gen_string(hand, YSTR(name)));
			}
		}
		YG(yajl_gen_array_close(hand));

		/* cohort names */
		YG(yajl_gen_string(hand, YSTR(COHORT_NAMES)));
		YG(yajl_gen_array_open(hand));
		{
			for (ic = 0; ic < nc; ic++) {
				const char *name = cohort_names[ic];
				YG(yajl_gen_string(hand, YSTR(name)));
			}
		}
		YG(yajl_gen_array_close(hand));

		/* coefficients */
		YG(yajl_gen_string(hand, YSTR(COEFFICIENTS)));
		const struct matrix *coefs = recv_fit_coefs(fit);
		YG(yajl_gen_matrix(hand, coefs));

		/* constraint names */
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

		/* constraints */
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
				YG(yajl_gen_integer(hand, dim * nc));

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

		/* constraint values */
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

		/* duals */
		YG(yajl_gen_string(hand, YSTR(DUALS)));
		const double *duals = recv_fit_duals(fit);
		YG(yajl_gen_vector(hand, ne, duals));

		/* count */
		YG(yajl_gen_string(hand, YSTR(COUNT)));
		YG(yajl_gen_array_open(hand));
		for (ic = 0; ic < nc; ic++) {
			size_t count = recv_loglik_count(ll, ic);
			YG(yajl_gen_integer(hand, count));
		}
		YG(yajl_gen_array_close(hand));

		/* score */
		YG(yajl_gen_string(hand, YSTR(SCORE)));
		YG(yajl_gen_array_open(hand));
		for (ic = 0; ic < nc; ic++) {
			const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
			const double *score = info->score;
			YG(yajl_gen_vector(hand, dim, score));
		}
		YG(yajl_gen_array_close(hand));

		/* information */
		YG(yajl_gen_string(hand, YSTR(INFORMATION)));
		YG(yajl_gen_array_open(hand));
		for (ic = 0; ic < nc; ic++) {
			const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
			const struct matrix *imat = &info->imat;
			YG(yajl_gen_matrix(hand, imat));
		}
		YG(yajl_gen_array_close(hand));

		/* rank */
		YG(yajl_gen_string(hand, YSTR(RANK)));
		size_t rank = dim * nc - ne;
		YG(yajl_gen_integer(hand, rank));

		/* deviance */
		YG(yajl_gen_string(hand, YSTR(DEVIANCE)));
		double dev = recv_fit_dev(fit);
		YG(yajl_gen_ieee754(hand, dev));

		/* df.residual */
		YG(yajl_gen_string(hand, YSTR(DF_RESIDUAL)));
		size_t ntot = recv_loglik_count_sum(ll);
		YG(yajl_gen_integer(hand, ntot - rank));

		/* null.deviance */
		YG(yajl_gen_string(hand, YSTR(NULL_DEVIANCE)));
		YG(yajl_gen_ieee754(hand, recv_fit_dev0(fit)));

		/* df.null */
		YG(yajl_gen_string(hand, YSTR(DF_NULL)));
		YG(yajl_gen_integer(hand, ntot));

		/* residuals */
		/* fitted.values */
		/* y = y */
		struct recv_resid resid;

		frame_clear(f);
		recv_resid_init(&resid, f, fit->ymsgs,
				fit->model.ncohort, fit->model.cohorts,
				coefs);
		YG(yajl_gen_string(hand, YSTR(OBSERVED_COUNTS)));
		YG(yajl_gen_matrix(hand, &resid.obs.dyad));

		YG(yajl_gen_string(hand, YSTR(EXPECTED_COUNTS)));
		YG(yajl_gen_matrix(hand, &resid.exp.dyad));

		recv_resid_deinit(&resid);
	}
	YG(yajl_gen_map_close(hand));

	/* effects */
	/* qr */
	/* linear.predictors (eta) */
	/* aic */
	/* iter */
	/* weights = wt */
	/* prior.weights = weights */
	/* converged = conv */
	/* boundary = boundary */
	/* call */
	/* formula */
	/* terms */
	/* data */
	/* offset */
	/* control */
	/* method */
	/* contrasts */
	/* xlevels */


	return err;
}


static void output(const struct recv_fit *fit)
{
	/* generator config */
	yajl_gen g = yajl_gen_alloc(NULL);
	yajl_gen_config(g, yajl_gen_beautify, 1);
	yajl_gen_recv_fit(g, fit);


	/* output */
	const unsigned char * buf;
	size_t len;
	yajl_gen_get_buf(g, &buf, &len);
	fwrite(buf, 1, len, stdout);
	yajl_gen_clear(g);
	yajl_gen_free(g);
}

static int do_fit(const struct messages *xmsgs, const struct messages *ymsgs,
		  const struct matrix *coefs0)
{
	size_t maxit = 100;
	size_t report = 1;
	bool trace = true;

	struct recv_fit fit;
	struct recv_fit_ctrl ctrl = RECV_FIT_CTRL0;
	ctrl.gtol = 1e-7;
	recv_fit_init(&fit, &frame, xmsgs, ymsgs, ncohort, cohorts, &ctrl);
	add_constraints(&fit);

	enum recv_fit_task task;
	size_t it = 0;

	//struct matrix coefs0;
	//matrix_init(&coefs0, dim, nc);
	//matrix_fill(&coefs0, 0.0001);

	for (it = 0, task = recv_fit_start(&fit, coefs0);
	     it < maxit && task == RECV_FIT_STEP;
	     it++, task = recv_fit_advance(&fit)) {
		if (trace && it % report == 0 && task != RECV_FIT_CONV) {
			// const struct recv_loglik *ll = recv_fit_loglik(&fit);
			// const struct recv_loglik_info *info = recv_loglik_info(ll);
			// size_t n = info->nrecv;
			// double dev = n * info->dev;
			// const struct vector *score = &info->score;
			// double ngrad = vector_max_abs(score);
			double step = recv_fit_step(&fit);
			double ngrad = recv_fit_grad_norm2(&fit);
			fprintf(stderr, "iter %"SSIZE_FMT"; |grad| = %.16f; step = %.16f\n",
				it, ngrad, step);

			// fprintf(stderr, "iter %"SSIZE_FMT" deviance = %.2f; |grad| = %.16f; step = %.16f\n",
			//	it, dev, ngrad, step);
		}
	}

	//matrix_deinit(&coefs0);

	int err = 0;

	if (task != RECV_FIT_CONV) {
		fprintf(stderr, "ERROR: %s\n", recv_fit_errmsg(&fit));
		err = -1;
	}

	output(&fit);
	recv_fit_deinit(&fit);

	return err;
}

static void init_coefs(struct matrix *coefs, char *filename)
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

	/* file read error handling */
	if (rd == 0 && !feof(fp)) {
		fprintf(stderr, "error encountered on file read\n");
		goto cleanup_file;
	}
	assert(rd == sz);
	filebuf[sz] = '\0';

	char errbuf[1024];
	errbuf[0] = '\0';

	yajl_val node = yajl_tree_parse(filebuf, errbuf, sizeof(errbuf));

	/* parse error handling */
	if (node == NULL) {
		fprintf(stderr, "parse_error: ");
		if (strlen(errbuf)) fprintf(stderr, " %s", errbuf);
		else fprintf(stderr, "unknown error");
		fprintf(stderr, "\n");
		goto cleanup_yajl;
	}

	/* ... and extract a nested value from the config file */
	const char * coefs_path[] = { COEFFICIENTS, NULL };
	const char * nrow_path[] = { "nrow", NULL};
	const char * ncol_path[] = { "ncol", NULL};
	const char * data_path[] = { "data", NULL};

	yajl_val vcoefs = yajl_tree_get(node, coefs_path, yajl_t_any);

	if (!vcoefs)
		goto cleanup_yajl;

	yajl_val vnrow = yajl_tree_get(vcoefs, nrow_path, yajl_t_number);
	yajl_val vncol = yajl_tree_get(vcoefs, ncol_path, yajl_t_number);
	yajl_val vdata = yajl_tree_get(vcoefs, data_path, yajl_t_array);

	if (!vnrow || !vncol || !vdata)
		goto cleanup_yajl;

	long long nrow = YAJL_GET_INTEGER(vnrow);
	long long ncol = YAJL_GET_INTEGER(vncol);
	size_t i, n = nrow * ncol;

	if (nrow < 0 || ncol < 0)
		goto cleanup_yajl;

	if (YAJL_GET_ARRAY(vdata)->len != n)
		goto cleanup_yajl;

	matrix_init(coefs, nrow, ncol);
	yajl_val *item = YAJL_GET_ARRAY(vdata)->values;
	double *ptr = matrix_to_ptr(coefs);

	for (i = 0; i < n; i++) {
		if (!YAJL_IS_NUMBER(item[i])) {
			matrix_deinit(coefs);
			err = 1;
			goto cleanup_coefs;
		}
		ptr[i] = YAJL_GET_DOUBLE(item[i]);
	}
	err = 0;
	goto cleanup_yajl;


cleanup_coefs:
	matrix_deinit(coefs);
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
	bool boot;
	int32_t seed;
};

static struct options parse_options(int argc, char **argv)
{
	int c;
	static struct option long_options[] = {
		{ "start", required_argument, 0, 's' },
		{ "boot", optional_argument, 0, 'b' },
		{ NULL, 0, NULL, 0 }
	};
	int option_index = 0;

	struct options opts;
	opts.startfile = NULL;
	opts.boot = false;
	opts.seed = 0;


	while ((c = getopt_long(argc, argv, "s:b:", long_options, &option_index)) != -1) {
		 switch (c) {
		 case 's':
			opts.startfile = optarg;
			fprintf(stderr, "start file: '%s'\n", opts.startfile);
			fflush(stderr);
			break;
		 case 'b':
			opts.boot = true;
			if (optarg) {
				opts.seed = (int32_t)strtol(optarg, NULL, 10);
			}
			fprintf(stderr, "bootstrap seed: '%" SSIZE_FMT"'\n", (size_t)opts.seed);
			fflush(stderr);
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
	struct matrix coefs0;
	bool has_coefs0 = 0;

	if (opts.startfile) {
		has_coefs0 = true;
		init_coefs(&coefs0, opts.startfile);
	}

	struct recv_boot boot;
	struct messages enron_messages;
	struct messages *xmsgs, *ymsgs;
	xmsgs = ymsgs = &enron_messages;

	setup();
	enron_messages_init(&enron_messages, 5);
	struct matrix *pcoefs0 = has_coefs0 ? &coefs0 : NULL;

	if (opts.boot) {
		dsfmt_t dsfmt;
		dsfmt_init_gen_rand(&dsfmt, opts.seed);

		recv_boot_init(&boot, &frame, &enron_messages, ncohort,
			       cohorts, pcoefs0, &dsfmt);
		ymsgs = &boot.messages;
	}

	err = do_fit(xmsgs, ymsgs, pcoefs0);

	messages_deinit(&enron_messages);

	if (opts.boot) {
		recv_boot_deinit(&boot);
	}

	teardown();

	if (has_coefs0) {
		matrix_deinit(&coefs0);
	}

	return err;
}
