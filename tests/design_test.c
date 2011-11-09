#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include "coreutil.h"
#include "cmockery.h"
#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"
#include "xalloc.h"

#include "enron.h"
#include "design.h"
#include "frame.h"


static size_t enron_nactor;
static size_t enron_ncohort;
static size_t enron_ntrait;
static size_t *enron_cohorts;
static struct dmatrix enron_traits;
static const char * const *enron_cohort_names;
static const char * const *enron_trait_names;

static struct frame frame;
static struct design *design;

static size_t nsend, nrecv;
static int has_loops;
static struct dmatrix traits;
static size_t ntrait;
static const char * const *trait_names;
static int has_effects;
static int has_loops;
static double *intvls;
size_t nintvl;
static size_t dim;


static double vector_dist(size_t n, double *x, double *y)
{
	size_t i;
	double dist = 0;

	for (i = 0; i < n; i++) {
		double d = x[i] - y[i];
		dist += d * d;
	}

	return sqrt(dist);
}

static void enron_setup_fixture()
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	enron_employees_init(&enron_nactor, &enron_cohorts, &enron_ncohort,
			     &enron_cohort_names, &enron_traits.data,
			     &enron_ntrait, &enron_trait_names);
	enron_traits.lda = MAX(1, enron_nactor);
}

static void enron_teardown_fixture()
{
	free(enron_cohorts);
	free(enron_traits.data);
	print_message("\n\n");
}

static void enron_setup()
{
	nsend = enron_nactor;
	nrecv = enron_nactor;
	has_loops = 0;
	traits = enron_traits;
	ntrait = enron_ntrait;
	trait_names = enron_trait_names;
	dim = enron_ntrait;
	has_effects = 0;
	has_loops = 0;
	intvls = NULL;
	nintvl = 0;

	frame_init(&frame, nsend, nrecv, has_loops, intvls, nintvl);
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
}

static void enron_teardown()
{
	frame_deinit(&frame);
	free(intvls);
}

static void enron_reff_setup_fixture()
{
	print_message("Enron employees (with receiver effects)\n");
	print_message("---------------------------------------\n");
	enron_employees_init(&enron_nactor, &enron_cohorts, &enron_ncohort,
			     &enron_cohort_names, &enron_traits.data,
			     &enron_ntrait, &enron_trait_names);
	enron_traits.lda = MAX(1, enron_nactor);
}

static void enron_reff_teardown_fixture()
{
	free(enron_cohorts);
	free(enron_traits.data);
	print_message("\n\n");
}

static void enron_reff_setup()
{
	nsend = enron_nactor;
	nrecv = enron_nactor;
	has_loops = 0;
	traits = enron_traits;
	ntrait = enron_ntrait;
	trait_names = enron_trait_names;
	dim = ntrait + nrecv;
	has_effects = 1;
	has_loops = 0;
	intvls = NULL;
	nintvl = 0;
	
	frame_init(&frame, nsend, nrecv, has_loops, intvls, nintvl);
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
}

static void enron_reff_teardown()
{
	frame_deinit(&frame);
	free(intvls);
}


static void test_size()
{
	assert_int_equal(design_dim(design), dim);
	assert_int_equal(design_count(design), nrecv);
}

static void matrix_assign_reffects(struct dmatrix *x,
				   size_t nrecv,
				   bool has_reffects)
{
	if (!has_reffects)
		return;

	size_t i;

	matrix_dzero(nrecv, nrecv, x);
	for (i = 0; i < nrecv; i++) {
		MATRIX_ITEM(x, i, i) = 1.0;
	}
}

static void matrix_init_design0(struct dmatrix *x, const struct design *d)
{
	const struct dmatrix *traits = design_traits(d);
	int has_reffects = design_has_effects(d);
	
	size_t nrecv = design_count(d);
	size_t pr = design_traits_dim(d);
	size_t ireff = 0;
	size_t nreff = has_reffects ? nrecv : 0;
	size_t istat = ireff + nreff;
	size_t nstat = pr;
	size_t dim =  istat + nstat;
	
	*x = (struct dmatrix) { xmalloc(nrecv * dim * sizeof(double)), MAX(1, nrecv) };
	struct dmatrix xreff = { MATRIX_COL(x, ireff), x->lda }; // nreff
	matrix_assign_reffects(&xreff, nrecv, has_reffects);
	struct dmatrix xstat = { MATRIX_COL(x, istat), x->lda };
	lapack_dlacpy(LA_COPY_ALL, nrecv, nstat, traits, &xstat);
}

static void test_mul0()
{

	struct dmatrix matrix;
	double *x, *y, *y1;	
	size_t i, n, p;
	
	n = design_count(design);
	p = design_dim(design);
	
	x = xmalloc(p * sizeof(double));
	y = xmalloc(n * sizeof(double));
	y1 = xmalloc(n * sizeof(double));
	
	matrix_init_design0(&matrix, design);
	
	for (i = 0; i < p; i++) {
		x[i] =  (7 * i) % 5 - 2.0;
	}
		
	design_mul0(1.0, BLAS_NOTRANS, design, x, 0.0, y);
	blas_dgemv(BLAS_NOTRANS, n, p, 1.0, &matrix, x, 1, 0.0, y1, 1);
	assert_true(vector_dist(n, y, y1) == 0.0);		
		
	design_mul0(1.0, BLAS_NOTRANS, design, x, 1.0, y);
	blas_dgemv(BLAS_NOTRANS, n, p, 1.0, &matrix, x, 1, 1.0, y1, 1);
	assert_true(vector_dist(n, y, y1) == 0.0);
		
	design_mul0(1.0, BLAS_NOTRANS, design, x, -1.0, y);
	blas_dgemv(BLAS_NOTRANS, n, p, 1.0, &matrix, x, 1, -1.0, y1, 1);	
	assert_true(vector_dist(n, y, y1) == 0.0);		
		
	design_mul0(2.0, BLAS_NOTRANS, design, x, 2.0, y);
	blas_dgemv(BLAS_NOTRANS, n, p, 2.0, &matrix, x, 1, 2.0, y1, 1);	
	assert_true(vector_dist(n, y, y1) == 0.0);		
		
	free(matrix.data);
	free(y1);
	free(y);
	free(x);
}


static void test_tmul0()
{
	
	struct dmatrix matrix;
	double *x, *y, *y1;	
	size_t i, n, p;
	
	n = design_count(design);
	p = design_dim(design);
	
	x = xmalloc(n * sizeof(double));
	y = xmalloc(p * sizeof(double));
	y1 = xmalloc(p * sizeof(double));
	matrix_init_design0(&matrix, design);
		
	for (i = 0; i < n; i++) {
		x[i] = (3 * i + 1) % 5 - 2.0;
	}
		
	design_mul0(1.0, BLAS_TRANS, design, x, 0.0, y);
	blas_dgemv(BLAS_TRANS, n, p, 1.0, &matrix, x, 1, 0.0, y1, 1);
	assert_true(vector_dist(p, y, y1) == 0.0);		
		
	design_mul0(1.0, BLAS_TRANS, design, x, 1.0, y);
	blas_dgemv(BLAS_TRANS, n, p, 1.0, &matrix, x, 1, 1.0, y1, 1);
	assert_true(vector_dist(p, y, y1) == 0.0);
		
	design_mul0(1.0, BLAS_TRANS, design, x, -1.0, y);
	blas_dgemv(BLAS_TRANS, n, p, 1.0, &matrix, x, 1, -1.0, y1, 1);	
	assert_true(vector_dist(p, y, y1) == 0.0);
		
	design_mul0(2.0, BLAS_TRANS, design, x, 2.0, y);
	blas_dgemv(BLAS_TRANS, n, p, 2.0, &matrix, x, 1, 2.0, y1, 1);
	assert_true(vector_dist(p, y, y1) == 0.0);
		
	free(matrix.data);
	free(y1);
	free(y);
	free(x);
}


static void test_tmuls0()
{
	struct dmatrix matrix;
	double *x, *y, *y1;
	struct vpattern pat;
	size_t i, n, p;

	n = design_count(design);
	p = design_dim(design);
	x = xmalloc(n * sizeof(x[0]));
	pat.indx = xmalloc(n * sizeof(pat.indx[0]));

	y = xmalloc(p * sizeof(double));
	y1 = xmalloc(p * sizeof(double));
	matrix_init_design0(&matrix, design);

	pat.nz = 0;
	for (i = 0; i < n; i++) {
		if (i % 3 == 0) {
			x[pat.nz] = (3 * i + 1) % 5 - 2.0;
			pat.indx[pat.nz] = i;
			pat.nz++;
		}
	}

	design_muls0(1.0, BLAS_TRANS, design, x, &pat, 0.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 1.0, &matrix, x, &pat, 0.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	design_muls0(1.0, BLAS_TRANS, design, x, &pat, 1.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 1.0, &matrix, x, &pat, 1.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	design_muls0(1.0, BLAS_TRANS, design, x, &pat, -1.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 1.0, &matrix, x, &pat, -1.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	design_muls0(2.0, BLAS_TRANS, design, x, &pat, 2.0, y);
	sblas_dgemvi(BLAS_TRANS, n, p, 2.0, &matrix, x, &pat, 2.0, y1);
	assert_true(vector_dist(p, y, y1) == 0.0);

	free(matrix.data);
	free(y1);
	free(y);
	free(pat.indx);
	free(x);
}


static void test_muls0()
{

	struct dmatrix matrix;
	double *x, *y, *y1;
	struct vpattern pat;
	size_t i, n, p;

	n = design_count(design);
	p = design_dim(design);
	x = xmalloc(n * sizeof(x[0]));
	pat.indx = xmalloc(n * sizeof(pat.indx[0]));

	y = xmalloc(n * sizeof(double));
	y1 = xmalloc(n * sizeof(double));
	matrix_init_design0(&matrix, design);

	pat.nz = 0;
	for (i = 0; i < p; i++) {
		if (i % 2 == 0) {
			x[pat.nz] = (7 * i) % 5 - 2.0;
			pat.indx[pat.nz] = i;
			pat.nz++;
		}
	}

	design_muls0(1.0, BLAS_NOTRANS, design, x, &pat, 0.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 1.0, &matrix, x, &pat, 0.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	design_muls0(1.0, BLAS_NOTRANS, design, x, &pat, 1.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 1.0, &matrix, x, &pat, 1.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	design_muls0(1.0, BLAS_NOTRANS, design, x, &pat, -1.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 1.0, &matrix, x, &pat, -1.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	design_muls0(2.0, BLAS_NOTRANS, design, x, &pat, 2.0, y);
	sblas_dgemvi(BLAS_NOTRANS, n, p, 2.0, &matrix, x, &pat, 2.0, y1);
	assert_true(vector_dist(n, y, y1) == 0.0);

	free(matrix.data);
	free(y1);
	free(y);
	free(pat.indx);
	free(x);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_size, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_mul0, enron_setup, enron_teardown),		
		unit_test_setup_teardown(test_tmul0, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_muls0, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_tmuls0, enron_setup, enron_teardown),	
		unit_test_teardown(enron_suite, enron_teardown_fixture),
		
		unit_test_setup(enron_reff_suite, enron_reff_setup_fixture),
		unit_test_setup_teardown(test_size, enron_reff_setup, enron_reff_teardown),
		unit_test_setup_teardown(test_mul0, enron_reff_setup, enron_reff_teardown),		
		unit_test_setup_teardown(test_tmul0, enron_reff_setup, enron_reff_teardown),
		unit_test_setup_teardown(test_muls0, enron_reff_setup, enron_reff_teardown),
		unit_test_setup_teardown(test_tmuls0, enron_reff_setup, enron_reff_teardown),	
		unit_test_teardown(enron_reff_suite, enron_reff_teardown_fixture),

	};
	return run_tests(tests);
}


