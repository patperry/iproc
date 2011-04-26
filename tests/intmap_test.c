#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "intmap.h"

static struct intmap map;
static intptr_t *keys;
static char *vals;
static ssize_t size;


static void empty_setup_fixture(void **state)
{
	print_message("empty intmap\n");
	print_message("------------\n");
}

static void teardown_fixture(void **state)
{
	print_message("\n\n");
}

static void empty_setup(void **state)
{
	static intptr_t *empty_keys = NULL;
	static char *empty_vals = NULL;

	intmap_init(&map, sizeof(char), alignof(char));
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

static void test_add(void **state)
{
	intptr_t key = 1;
	char val = 'a';
	intmap_add(&map, key, &val);
	assert_int_equal(intmap_size(&map), size + 1);
	assert_true(intmap_contains(&map, key));
	assert_int_equal(*(char *)intmap_lookup(&map, key), val);
}

static void test_add_hard(void **state)
{
	char val = 'a';
	int i, j, n = 100;
	
	for (i = 0; i < n; i++) {
		val = (char)(5 * i + 1);
		intmap_add(&map, i, &val);
		assert_int_equal(intmap_size(&map), size + i + 1);
		assert_true(intmap_contains(&map, i));
		assert_int_equal(*(char *)intmap_lookup(&map, i), val);
		
		for (j = 0; j < i; j++) {
			assert_true(intmap_contains(&map, j));
			assert_int_equal(*(char *)intmap_lookup(&map, j),
					 (char)(5 * j + 1));
			
		}
	}
}


int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_size, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_add, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_add_hard, empty_setup, empty_teardown),		
		unit_test_teardown(empty_suite, teardown_fixture),
	};
	return run_tests(tests);
}


