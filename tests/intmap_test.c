#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "intmap.h"

static struct intmap map;
static intptr_t *keys;
static int *vals;
static ssize_t size;


static void empty_setup_fixture(void **state)
{
	print_message("empty intset\n");
	print_message("------------\n");
}

static void teardown_fixture(void **state)
{
	print_message("\n\n");
}

static void empty_setup(void **state)
{
	static intptr_t empty_keys[] = { };
	static int empty_vals[] = { };

	intmap_init(&map, int);
	keys = empty_keys;
	vals = empty_vals;
	size = 0;
}

static void empty_teardown(void **state)
{
	intmap_deinit(&map);
}

static void test_size(void **state)
{
	assert_int_equal(intmap_size(&map), size);
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_size, empty_setup, empty_teardown),
		unit_test_teardown(empty_suite, teardown_fixture),
	};
	return run_tests(tests);
}


