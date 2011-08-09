#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmockery.h>

#include "enron.h"
#include "actors.h"

static struct actors actors;
static struct matrix traits;
static struct matrix matrix;

static void enron_setup_fixture(void **state)
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	enron_employees_init(&actors, &traits);
	
	ssize_t n = actors_count(&actors);
	ssize_t c, nc = actors_cohort_count(&actors);
	matrix_init(&matrix, n, nc);
	
	const struct cohort *cohort = actors_cohorts(&actors);
	
	for (c = 0; c < nc; c++) {
		const ssize_t *items = array_to_ptr(&cohort[c].actors);
		ssize_t ia, na = array_count(&cohort[c].actors);
		for (ia = 0; ia < na; ia++) {
			ssize_t i = items[ia];
			matrix_set_item(&matrix, i, c, 1.0);
		}
	}
}

static void enron_teardown_fixture(void **state)
{
	matrix_deinit(&matrix);
	matrix_deinit(&traits);
	actors_deinit(&actors);	
	print_message("\n\n");
}

static void enron_setup(void **state)
{
}

static void enron_teardown(void **state)
{
}

static void enron_test_size(void **state)
{
	assert_int_equal(actors_count(&actors), matrix_nrow(&matrix));
	assert_int_equal(actors_cohort_count(&actors), matrix_nrow(&traits));	
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
		vector_set_item(&x, i, (7 * i) % 5 - 2);
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
		vector_set_item(&x, i, (3 * i + 1) % 5 - 2);
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
	struct svector x, y;
	struct vector y0, diff, zero;
	ssize_t n = matrix_nrow(&matrix);
	ssize_t p = matrix_ncol(&matrix);
	ssize_t i;
	
	svector_init(&x, p);
	svector_init(&y, n);
	vector_init(&y0, n);
	vector_init(&diff, n);
	vector_init(&zero, n);	
	vector_fill(&zero, 0.0);
	
	for (i = 0; i < p; i++) {
		if (i % 2 == 0)
			svector_set_item(&x, i, (7 * i) % 5 - 2);
	}
	
	actors_muls(1.0, TRANS_NOTRANS, &actors, &x, 0.0, &y);
	matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, 0.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	actors_muls(1.0, TRANS_NOTRANS, &actors, &x, 1.0, &y);
	matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, 1.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	actors_muls(1.0, TRANS_NOTRANS, &actors, &x, -1.0, &y);
	matrix_muls(1.0, TRANS_NOTRANS, &matrix, &x, -1.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	actors_muls(2.0, TRANS_NOTRANS, &actors, &x, 2.0, &y);
	matrix_muls(2.0, TRANS_NOTRANS, &matrix, &x, 2.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	vector_deinit(&zero);
	vector_deinit(&diff);
	vector_deinit(&y0);
	svector_deinit(&y);
	svector_deinit(&x);
}

static void test_tmuls(void **state)
{
	struct svector x, y;
	struct vector y0, diff, zero;
	ssize_t n = matrix_nrow(&matrix);
	ssize_t p = matrix_ncol(&matrix);
	ssize_t i;
	
	svector_init(&x, n);
	svector_init(&y, p);
	vector_init(&y0, p);
	vector_init(&diff, p);	
	vector_init(&zero, p);	
	
	for (i = 0; i < n; i++) {
		if (i % 3 == 0)
			svector_set_item(&x, i, (3 * i + 1) % 5 - 2);
	}
	
	actors_muls(1.0, TRANS_TRANS, &actors, &x, 0.0, &y);
	matrix_muls(1.0, TRANS_TRANS, &matrix, &x, 0.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	actors_muls(1.0, TRANS_TRANS, &actors, &x, 1.0, &y);
	matrix_muls(1.0, TRANS_TRANS, &matrix, &x, 1.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	actors_muls(1.0, TRANS_TRANS, &actors, &x, -1.0, &y);
	matrix_muls(1.0, TRANS_TRANS, &matrix, &x, -1.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	actors_muls(2.0, TRANS_TRANS, &actors, &x, 2.0, &y);
	matrix_muls(2.0, TRANS_TRANS, &matrix, &x, 2.0, &y0);
	vector_assign_copy(&diff, &y0);
	svector_axpy(-1.0, &y, &diff);
	assert_true(vector_equals(&diff, &zero));
	
	vector_deinit(&zero);
	vector_deinit(&diff);
	vector_deinit(&y0);
	svector_deinit(&y);
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


