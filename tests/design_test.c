#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include "cmockery.h"

#include "enron.h"
#include "design.h"


static size_t enron_nactor;
static size_t enron_ncohort;
static ptrdiff_t *enron_cohorts;
static struct matrix enron_traits;
static const char * const *enron_cohort_names;
static const char * const *enron_trait_names;

static struct design design;

static size_t nsend, nrecv;
static struct matrix recv_traits;
static bool has_reffects;
static bool has_loops;
static struct vector intervals;
static ssize_t dim;


static void enron_setup_fixture()
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	enron_employees_init(&enron_nactor, &enron_ncohort, &enron_cohorts,
			     &enron_cohort_names, &enron_traits,
			     &enron_trait_names);
}

static void enron_teardown_fixture()
{
	free(enron_cohorts);
	matrix_deinit(&enron_traits);
	print_message("\n\n");
}

static void enron_setup()
{
	nsend = enron_nactor;
	nrecv = enron_nactor;
	matrix_init_copy(&recv_traits, BLAS_NOTRANS, &enron_traits);
	dim = matrix_ncol(&recv_traits);
	has_reffects = false;
	has_loops = false;
	vector_init(&intervals, 0);

	design_init(&design, nsend, nrecv, &recv_traits, enron_trait_names, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
}

static void enron_teardown()
{
	design_deinit(&design);
	vector_deinit(&intervals);
	matrix_deinit(&recv_traits);
}

static void enron_reff_setup_fixture()
{
	print_message("Enron employees (with receiver effects)\n");
	print_message("---------------------------------------\n");
	enron_employees_init(&enron_nactor, &enron_ncohort, &enron_cohorts,
			     &enron_cohort_names, &enron_traits,
			     &enron_trait_names);
}

static void enron_reff_teardown_fixture()
{
	free(enron_cohorts);
	matrix_deinit(&enron_traits);
	print_message("\n\n");
}

static void enron_reff_setup()
{
	nsend = enron_nactor;
	nrecv = enron_nactor;
	matrix_init_copy(&recv_traits, BLAS_NOTRANS, &enron_traits);
	dim = matrix_ncol(&recv_traits) + nrecv;
	has_reffects = true;
	has_loops = false;
	vector_init(&intervals, 0);
	
	design_init(&design, nsend, nrecv, &recv_traits, enron_trait_names, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
}

static void enron_reff_teardown()
{
	design_deinit(&design);
	vector_deinit(&intervals);
	matrix_deinit(&recv_traits);
}


static void test_size()
{
	assert_int_equal(design_recv_dim(&design), dim);
	assert_int_equal(design_send_count(&design), nsend);
	assert_int_equal(design_recv_count(&design), nrecv);
}

static void matrix_assign_static(struct matrix *x,
				 const struct matrix *traits)
{
	matrix_assign_copy(x, BLAS_NOTRANS, traits);
}

static void matrix_assign_reffects(struct matrix *x,
				   size_t nrecv,
				   bool has_reffects)
{
	(void)nrecv;
	assert((size_t)matrix_nrow(x) == nrecv);

	if (has_reffects) {
		assert((size_t)matrix_ncol(x) == nrecv);
		matrix_assign_identity(x);
	} else {
		assert(matrix_ncol(x) == 0);
	}
}

static void matrix_init_design0(struct matrix *x, const struct design *d)
{
	const struct matrix *traits = design_traits(d);
	bool has_reffects = design_recv_effects(d);
	struct matrix xstat, xreff;
	
	ssize_t nrecv = design_recv_count(d);
	ssize_t pr = matrix_ncol(traits);
	ssize_t ireff = 0;
	ssize_t nreff = has_reffects ? nrecv : 0;
	ssize_t istat = ireff + nreff;
	ssize_t nstat = pr;
	ssize_t dim =  istat + nstat;
	
	matrix_init(x, nrecv, dim);
	xreff = matrix_slice_cols(x, ireff, nreff);
	matrix_assign_reffects(&xreff, nrecv, has_reffects);
	xstat = matrix_slice_cols(x, istat, nstat);
	matrix_assign_static(&xstat, traits);
}

static void test_mul0()
{

	struct matrix matrix;
	struct vector x, y, y1;	
	ssize_t i, n, p;
	
	n = design_recv_count(&design);
	p = design_recv_dim(&design);
	
	vector_init(&x, p);
	vector_init(&y, n);
	vector_init(&y1, n);
	
	matrix_init_design0(&matrix, &design);
	
	for (i = 0; i < p; i++) {
		vector_set_item(&x, i, (7 * i) % 5 - 2);
	}
		
	design_recv_mul0(1.0, BLAS_NOTRANS, &design, &x, 0.0, &y);
	matrix_mul(1.0, BLAS_NOTRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	design_recv_mul0(1.0, BLAS_NOTRANS, &design, &x, 1.0, &y);
	matrix_mul(1.0, BLAS_NOTRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	design_recv_mul0(1.0, BLAS_NOTRANS, &design, &x, -1.0, &y);
	matrix_mul(1.0, BLAS_NOTRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	design_recv_mul0(2.0, BLAS_NOTRANS, &design, &x, 2.0, &y);
	matrix_mul(2.0, BLAS_NOTRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	matrix_deinit(&matrix);
	
	vector_deinit(&y1);
	vector_deinit(&y);
	vector_deinit(&x);
}


static void test_tmul0()
{
	
	struct matrix matrix;
	struct vector x, y, y1;	
	ssize_t i, n, p;
	
	n = design_recv_count(&design);
	p = design_recv_dim(&design);
	
	vector_init(&x, n);
	vector_init(&y, p);
	vector_init(&y1, p);
	matrix_init_design0(&matrix, &design);
		
	for (i = 0; i < n; i++) {
		vector_set_item(&x, i, (3 * i + 1) % 5 - 2);
	}
		
	design_recv_mul0(1.0, BLAS_TRANS, &design, &x, 0.0, &y);
	matrix_mul(1.0, BLAS_TRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	design_recv_mul0(1.0, BLAS_TRANS, &design, &x, 1.0, &y);
	matrix_mul(1.0, BLAS_TRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	design_recv_mul0(1.0, BLAS_TRANS, &design, &x, -1.0, &y);
	matrix_mul(1.0, BLAS_TRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	design_recv_mul0(2.0, BLAS_TRANS, &design, &x, 2.0, &y);
	matrix_mul(2.0, BLAS_TRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	matrix_deinit(&matrix);
	vector_deinit(&y1);
	vector_deinit(&y);
	vector_deinit(&x);
}


static void test_tmuls0()
{
	
	struct matrix matrix;
	struct svector x;
	struct vector y, y1;	
	ssize_t i, n, p;
	
	n = design_recv_count(&design);
	p = design_recv_dim(&design);
	
	svector_init(&x, n);
	vector_init(&y, p);
	vector_init(&y1, p);
	matrix_init_design0(&matrix, &design);
		
	for (i = 0; i < n; i++) {
		if (i % 3 == 0)
			svector_set_item(&x, i, (3 * i + 1) % 5 - 2);
	}
		
	design_recv_muls0(1.0, BLAS_TRANS, &design, &x, 0.0, &y);
	matrix_muls(1.0, BLAS_TRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	design_recv_muls0(1.0, BLAS_TRANS, &design, &x, 1.0, &y);
	matrix_muls(1.0, BLAS_TRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	design_recv_muls0(1.0, BLAS_TRANS, &design, &x, -1.0, &y);
	matrix_muls(1.0, BLAS_TRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	design_recv_muls0(2.0, BLAS_TRANS, &design, &x, 2.0, &y);
	matrix_muls(2.0, BLAS_TRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	matrix_deinit(&matrix);
	vector_deinit(&y1);
	vector_deinit(&y);
	svector_deinit(&x);
}


static void test_muls0()
{
	
	struct matrix matrix;
	struct svector x;
	struct vector y, y1;	
	ssize_t i, n, p;
	
	n = design_recv_count(&design);
	p = design_recv_dim(&design);
	
	svector_init(&x, p);
	vector_init(&y, n);
	vector_init(&y1, n);
	matrix_init_design0(&matrix, &design);

	for (i = 0; i < p; i++) {
		if (i % 2 == 0)
			svector_set_item(&x, i, (7 * i) % 5 - 2);
	}
		
	design_recv_muls0(1.0, BLAS_NOTRANS, &design, &x, 0.0, &y);
	matrix_muls(1.0, BLAS_NOTRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	design_recv_muls0(1.0, BLAS_NOTRANS, &design, &x, 1.0, &y);
	matrix_muls(1.0, BLAS_NOTRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);
		
	design_recv_muls0(1.0, BLAS_NOTRANS, &design, &x, -1.0, &y);
	matrix_muls(1.0, BLAS_NOTRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	design_recv_muls0(2.0, BLAS_NOTRANS, &design, &x, 2.0, &y);
	matrix_muls(2.0, BLAS_NOTRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_dist(&y, &y1) == 0.0);		
		
	matrix_deinit(&matrix);
	vector_deinit(&y1);
	vector_deinit(&y);
	svector_deinit(&x);
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


