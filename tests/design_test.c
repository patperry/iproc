#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include "cmockery.h"

#include "enron.h"
#include "design.h"


static struct matrix enron_employees;

static struct design design;

static struct actors senders;
static struct actors receivers;
static bool has_reffects;
static ssize_t dim;


static void enron_setup_fixture(void **state)
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	enron_employee_matrix_init(&enron_employees);
}

static void enron_teardown_fixture(void **state)
{
	matrix_deinit(&enron_employees);
	print_message("\n\n");
}

static void enron_setup(void **state)
{
	struct matrix enron_employees0;

	matrix_init_slice_cols(&enron_employees0, &enron_employees, 1,
			       matrix_ncol(&enron_employees) - 1);
	actors_init_matrix(&senders, &enron_employees, TRANS_NOTRANS);
	actors_init_matrix(&receivers, &enron_employees0, TRANS_NOTRANS);	
	dim = actors_dim(&senders) * actors_dim(&receivers);
	has_reffects = false;

	design_init(&design, &senders, &receivers, has_reffects);
	matrix_deinit(&enron_employees0);
}

static void enron_teardown(void **state)
{
	actors_deinit(&receivers);
	actors_deinit(&senders);
}

static void enron_reff_setup_fixture(void **state)
{
	print_message("Enron employees (with receiver effects)\n");
	print_message("---------------------------------------\n");
	enron_employee_matrix_init(&enron_employees);
}

static void enron_reff_teardown_fixture(void **state)
{
	matrix_deinit(&enron_employees);
	print_message("\n\n");
}

static void enron_reff_setup(void **state)
{
	struct matrix enron_employees0;
	
	matrix_init_slice_cols(&enron_employees0, &enron_employees, 1,
			       matrix_ncol(&enron_employees) - 1);
	actors_init_matrix(&senders, &enron_employees, TRANS_NOTRANS);
	actors_init_matrix(&receivers, &enron_employees0, TRANS_NOTRANS);	
	dim = (actors_size(&receivers)
	       + actors_dim(&senders) * actors_dim(&receivers));
	has_reffects = true;
	
	design_init(&design, &senders, &receivers, has_reffects);
	matrix_deinit(&enron_employees0);
}

static void enron_reff_teardown(void **state)
{
	actors_deinit(&receivers);
	actors_deinit(&senders);
}


static void test_size(void **state)
{
	assert_int_equal(design_dim(&design), dim);
	assert_int_equal(design_nsender(&design), actors_size(&senders));
	assert_int_equal(design_nreceiver(&design), actors_size(&receivers));
}

static void matrix_assign_static(struct matrix *x,
				 const struct actors *s,
				 const struct actors *r,
				 ssize_t isend)
{
	ssize_t irecv, jvar, nrecv, js, jr, ps, pr;
	double sj, rj;
	
	nrecv = actors_size(r);
	ps = actors_dim(s);
	pr = actors_dim(r);
	
	assert(matrix_nrow(x) == nrecv);
	assert(matrix_ncol(x) == ps * pr);

	for (js = 0; js < ps; js++) {
		sj = actors_get(s, isend, js);
		for (jr = 0; jr < pr; jr++) {
			jvar = js + jr * ps;			
			for (irecv = 0; irecv < nrecv; irecv++) {
				rj = actors_get(r, irecv, jr);
				matrix_set(x, irecv, jvar, sj * rj);
			}
		}
	}
}

static void matrix_assign_reffects(struct matrix *x,
				   const struct actors *r,
				   bool has_reffects)
{
	assert(matrix_nrow(x) == actors_size(r));	

	if (has_reffects) {
		assert(matrix_ncol(x) == actors_size(r));
		matrix_assign_identity(x);
	} else {
		assert(matrix_ncol(x) == 0);
	}
}

static void matrix_init_design0(struct matrix *x, const struct design *d,
				ssize_t isend)
{
	assert(0 <= isend && isend < design_nsender(d));
	
	const struct actors *s = design_senders(d);
	const struct actors *r = design_receivers(d);
	bool has_reffects = design_has_reffects(d);
	struct matrix xstat, xreff;
	
	ssize_t nrecv = actors_size(r);
	ssize_t ps = actors_dim(s);
	ssize_t pr = actors_dim(r);
	ssize_t ireff = 0;
	ssize_t nreff = has_reffects ? nrecv : 0;
	ssize_t istat = ireff + nreff;
	ssize_t nstat = ps * pr;
	ssize_t dim =  istat + nstat;
	
	matrix_init(x, nrecv, dim);
	matrix_init_slice_cols(&xreff, x, ireff, nreff);
	matrix_assign_reffects(&xreff, r, has_reffects);
	matrix_init_slice_cols(&xstat, x, istat, nstat);
	matrix_assign_static(&xstat, s, r, isend);
}

static void test_mul0(void **state)
{

	struct matrix matrix;
	struct vector x, y, y1;	
	ssize_t isend, nsend, i, n, p;
	
	n = design_nreceiver(&design);
	p = design_dim(&design);
	
	vector_init(&x, p);
	vector_init(&y, n);
	vector_init(&y1, n);
	
	nsend = actors_size(&senders);
	for (isend = 0; isend < nsend; isend += 10) {
		matrix_init_design0(&matrix, &design, isend);
	
		for (i = 0; i < p; i++) {
			vector_set(&x, i, (7 * i) % 5 - 2);
		}
		
		design_mul0(1.0, TRANS_NOTRANS, &design, isend, &x, 0.0, &y);
		matrix_mul(1.0, TRANS_NOTRANS, &matrix, &x, 0.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		
		design_mul0(1.0, TRANS_NOTRANS, &design, isend, &x, 1.0, &y);
		matrix_mul(1.0, TRANS_NOTRANS, &matrix, &x, 1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		design_mul0(1.0, TRANS_NOTRANS, &design, isend, &x, -1.0, &y);
		matrix_mul(1.0, TRANS_NOTRANS, &matrix, &x, -1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		assert_true(vector_equals(&y, &y1));
		
		design_mul0(2.0, TRANS_NOTRANS, &design, isend, &x, 2.0, &y);
		matrix_mul(2.0, TRANS_NOTRANS, &matrix, &x, 2.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		
		matrix_deinit(&matrix);
	}
	
	vector_deinit(&y1);
	vector_deinit(&y);
	vector_deinit(&x);
}


static void test_tmul0(void **state)
{
	
	struct matrix matrix;
	struct vector x, y, y1;	
	ssize_t isend, nsend, i, n, p;
	
	n = design_nreceiver(&design);
	p = design_dim(&design);
	
	vector_init(&x, n);
	vector_init(&y, p);
	vector_init(&y1, p);
	
	nsend = actors_size(&senders);
	for (isend = 1; isend < nsend; isend += 11) {
		matrix_init_design0(&matrix, &design, isend);
		
		for (i = 0; i < n; i++) {
			vector_set(&x, i, (3 * i + 1) % 5 - 2);
		}
		
		design_mul0(1.0, TRANS_TRANS, &design, isend, &x, 0.0, &y);
		matrix_mul(1.0, TRANS_TRANS, &matrix, &x, 0.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		
		design_mul0(1.0, TRANS_TRANS, &design, isend, &x, 1.0, &y);
		matrix_mul(1.0, TRANS_TRANS, &matrix, &x, 1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		design_mul0(1.0, TRANS_TRANS, &design, isend, &x, -1.0, &y);
		matrix_mul(1.0, TRANS_TRANS, &matrix, &x, -1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		design_mul0(2.0, TRANS_TRANS, &design, isend, &x, 2.0, &y);
		matrix_mul(2.0, TRANS_TRANS, &matrix, &x, 2.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		matrix_deinit(&matrix);
	}
	
	vector_deinit(&y1);
	vector_deinit(&y);
	vector_deinit(&x);
}


static void test_tmuls0(void **state)
{
	
	struct matrix matrix;
	struct svector x;
	struct vector y, y1;	
	ssize_t isend, nsend, i, n, p;
	
	n = design_nreceiver(&design);
	p = design_dim(&design);
	
	svector_init(&x, n);
	vector_init(&y, p);
	vector_init(&y1, p);
	
	nsend = actors_size(&senders);
	for (isend = 2; isend < nsend; isend += 12) {
		matrix_init_design0(&matrix, &design, isend);
		
		for (i = 0; i < n; i++) {
			if (i % 3 == 0)
				svector_set(&x, i, (3 * i + 1) % 5 - 2);
		}
		
		design_muls0(1.0, TRANS_TRANS, &design, isend, &x, 0.0, &y);
		matrix_muls(1.0, TRANS_TRANS, &matrix, &x, 0.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		
		design_muls0(1.0, TRANS_TRANS, &design, isend, &x, 1.0, &y);
		matrix_muls(1.0, TRANS_TRANS, &matrix, &x, 1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		design_muls0(1.0, TRANS_TRANS, &design, isend, &x, -1.0, &y);
		matrix_muls(1.0, TRANS_TRANS, &matrix, &x, -1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		design_muls0(2.0, TRANS_TRANS, &design, isend, &x, 2.0, &y);
		matrix_muls(2.0, TRANS_TRANS, &matrix, &x, 2.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		matrix_deinit(&matrix);
	}
	
	vector_deinit(&y1);
	vector_deinit(&y);
	svector_deinit(&x);
}


static void test_muls0(void **state)
{
	
	struct matrix matrix;
	struct svector x;
	struct vector y, y1;	
	ssize_t isend, nsend, i, n, p;
	
	n = design_nreceiver(&design);
	p = design_dim(&design);
	
	svector_init(&x, p);
	vector_init(&y, n);
	vector_init(&y1, n);
	
	nsend = actors_size(&senders);
	for (isend = 3; isend < nsend; isend += 13) {
		matrix_init_design0(&matrix, &design, isend);

		for (i = 0; i < p; i++) {
			if (i % 2 == 0)
				svector_set(&x, i, (7 * i) % 5 - 2);
		}
		
		design_muls0(1.0, TRANS_NOTRANS, &design, isend, &x, 0.0, &y);
		matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, 0.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		
		design_muls0(1.0, TRANS_NOTRANS, &design, isend, &x, 1.0, &y);
		matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, 1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);
		
		design_muls0(1.0, TRANS_NOTRANS, &design, isend, &x, -1.0, &y);
		matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, -1.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		assert_true(vector_equals(&y, &y1));
		
		design_muls0(2.0, TRANS_NOTRANS, &design, isend, &x, 2.0, &y);
		matrix_muls(2.0, TRANS_NOTRANS, &matrix, &x, 2.0, &y1);
		assert_true(vector_dist(&y, &y1) == 0.0);		
		
		matrix_deinit(&matrix);
	}
	
	vector_deinit(&y1);
	vector_deinit(&y);
	svector_deinit(&x);
}



int main(int argc, char **argv)
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


