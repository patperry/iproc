#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "array.h"

static struct array array;
static ssize_t size;
static ssize_t *elts;		// elements sorted in decending order

static void empty_setup_fixture()
{
	print_message("empty array\n");
	print_message("-----------\n");
}

static void empty_setup()
{
	static ssize_t *empty_elts = NULL;

	array_init(&array, sizeof(ssize_t));
	size = 0;
	elts = empty_elts;
}

static void singleton_setup_fixture()
{
	print_message("singleton array\n");
	print_message("---------------\n");
}

static void singleton_setup()
{
	static ssize_t singleton_elts[] = { 1234 };

	array_init(&array, sizeof(ssize_t));
	elts = singleton_elts;
	size = 1;
	array_add_range(&array, singleton_elts, size);
}

static void teardown()
{
	array_deinit(&array);
}

static void teardown_fixture()
{
	print_message("\n\n");
}

static void test_size()
{
	assert_int_equal(array_count(&array), size);
}

static void test_insert()
{
	ssize_t i, j;
	ssize_t val = 31337;

	for (i = 0; i <= size; i++) {
		struct array a;

		array_init_copy(&a, &array);

		array_insert(&a, i, &val);

		assert_int_equal(array_count(&a), size + 1);
		for (j = 0; j <= size; j++) {
			if (j < i) {
				assert_int_equal(*(ssize_t *)array_item(&a, j),
						 elts[j]);
			} else if (j == i) {
				assert_int_equal(*(ssize_t *)array_item(&a, j),
						 val);
			} else {
				assert_int_equal(*(ssize_t *)array_item(&a, j),
						 elts[j - 1]);
			}
		}

		array_deinit(&a);
	}
}

int main()
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_size, empty_setup, teardown),
		unit_test_setup_teardown(test_insert, empty_setup, teardown),
		unit_test_teardown(empty_suite, teardown_fixture),

		unit_test_setup(singleton_suite, singleton_setup_fixture),
		unit_test_setup_teardown(test_size, singleton_setup, teardown),
		unit_test_setup_teardown(test_insert, singleton_setup, teardown),
		unit_test_teardown(singleton_suite, teardown_fixture),
	};
	return run_tests(tests);
}
