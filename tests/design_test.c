#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <cmockery.h>

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
	actors_init_matrix(&senders, &enron_employees0, TRANS_NOTRANS);
	actors_init_matrix(&receivers, &enron_employees, TRANS_NOTRANS);
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

static void test_size(void **state)
{
	assert_int_equal(design_dim(&design), dim);
	assert_int_equal(design_nsender(&design), actors_size(&senders));
	assert_int_equal(design_nreceiver(&design), actors_size(&receivers));
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_size, enron_setup, enron_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}


