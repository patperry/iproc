#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmockery.h>

#include "enron.h"
#include "actors.h"


static struct matrix employees;
static struct matrix matrix;
static struct actors actors;


static void enron_setup_fixture(void **state)
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	enron_employee_matrix_init(&employees);
}

static void enron_teardown_fixture(void **state)
{
	matrix_deinit(&employees);
	print_message("\n\n");
}

static void enron_setup(void **state)
{
	matrix_init_copy(&matrix, &employees);
	actors_init_matrix(&actors, &matrix, TRANS_NOTRANS);
}

static void enron_teardown(void **state)
{
	matrix_deinit(&matrix);
	actors_deinit(&actors);
}

static void enron_test_size(void **state)
{
	assert_int_equal(actors_size(&actors), matrix_nrow(&matrix));
	assert_int_equal(actors_dim(&actors), matrix_ncol(&matrix));
	assert_int_equal(actors_cohorts_size(&actors), 12);
}

static void test_mul(void **state)
{
	struct vector x, y, y1;
	ssize_t n = matrix_nrow(&matrix);
	ssize_t p = matrix_ncol(&matrix);
	ssize_t i;
	
	vector_init(&x, p);
	vector_init(&y, n);
	vector_init(&y1, n);
	
	for (i = 0; i < p; i++) {
		vector_set(&x, i, (7 * i) % 5 - 2);
	}
	
	actors_mul(1.0, TRANS_NOTRANS, &actors, &x, 0.0, &y);
	matrix_mul(1.0, TRANS_NOTRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_equals(&y, &y1));

	actors_mul(1.0, TRANS_NOTRANS, &actors, &x, 1.0, &y);
	matrix_mul(1.0, TRANS_NOTRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_mul(1.0, TRANS_NOTRANS, &actors, &x, -1.0, &y);
	matrix_mul(1.0, TRANS_NOTRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_equals(&y, &y1));

	actors_mul(2.0, TRANS_NOTRANS, &actors, &x, 2.0, &y);
	matrix_mul(2.0, TRANS_NOTRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	vector_deinit(&y1);
	vector_deinit(&y);
	vector_deinit(&x);
}

static void test_tmul(void **state)
{
	struct vector x, y, y1;
	ssize_t n = matrix_nrow(&matrix);
	ssize_t p = matrix_ncol(&matrix);
	ssize_t i;
	
	vector_init(&x, n);
	vector_init(&y, p);
	vector_init(&y1, p);
	
	for (i = 0; i < n; i++) {
		vector_set(&x, i, (3 * i + 1) % 5 - 2);
	}
	
	actors_mul(1.0, TRANS_TRANS, &actors, &x, 0.0, &y);
	matrix_mul(1.0, TRANS_TRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_mul(1.0, TRANS_TRANS, &actors, &x, 1.0, &y);
	matrix_mul(1.0, TRANS_TRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_mul(1.0, TRANS_TRANS, &actors, &x, -1.0, &y);
	matrix_mul(1.0, TRANS_TRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_mul(2.0, TRANS_TRANS, &actors, &x, 2.0, &y);
	matrix_mul(2.0, TRANS_TRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	vector_deinit(&y1);
	vector_deinit(&y);
	vector_deinit(&x);
}

static void test_muls(void **state)
{
	struct svector x;
	struct vector y, y1;
	ssize_t n = matrix_nrow(&matrix);
	ssize_t p = matrix_ncol(&matrix);
	ssize_t i;
	
	svector_init(&x, p);
	vector_init(&y, n);
	vector_init(&y1, n);
	
	for (i = 0; i < p; i++) {
		if (i % 2 == 0)
			svector_set(&x, i, (7 * i) % 5 - 2);
	}
	
	actors_muls(1.0, TRANS_NOTRANS, &actors, &x, 0.0, &y);
	matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_muls(1.0, TRANS_NOTRANS, &actors, &x, 1.0, &y);
	matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_muls(1.0, TRANS_NOTRANS, &actors, &x, -1.0, &y);
	matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_muls(2.0, TRANS_NOTRANS, &actors, &x, 2.0, &y);
	matrix_muls(2.0, TRANS_NOTRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	vector_deinit(&y1);
	vector_deinit(&y);
	svector_deinit(&x);
}

static void test_tmuls(void **state)
{
	struct svector x;
	struct vector y, y1;
	ssize_t n = matrix_nrow(&matrix);
	ssize_t p = matrix_ncol(&matrix);
	ssize_t i;
	
	svector_init(&x, n);
	vector_init(&y, p);
	vector_init(&y1, p);
	
	for (i = 0; i < n; i++) {
		if (i % 3 == 0)
			svector_set(&x, i, (3 * i + 1) % 5 - 2);
	}
	
	actors_muls(1.0, TRANS_TRANS, &actors, &x, 0.0, &y);
	matrix_muls(1.0, TRANS_TRANS, &matrix, &x, 0.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_muls(1.0, TRANS_TRANS, &actors, &x, 1.0, &y);
	matrix_muls(1.0, TRANS_TRANS, &matrix, &x, 1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_muls(1.0, TRANS_TRANS, &actors, &x, -1.0, &y);
	matrix_muls(1.0, TRANS_TRANS, &matrix, &x, -1.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	actors_muls(2.0, TRANS_TRANS, &actors, &x, 2.0, &y);
	matrix_muls(2.0, TRANS_TRANS, &matrix, &x, 2.0, &y1);
	assert_true(vector_equals(&y, &y1));
	
	vector_deinit(&y1);
	vector_deinit(&y);
	svector_deinit(&x);
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(enron_test_size, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_mul, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_tmul, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_muls, enron_setup, enron_teardown),		
		unit_test_setup_teardown(test_tmuls, enron_setup, enron_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}


