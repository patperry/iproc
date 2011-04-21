#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "hashset.h"

DEFINE_HASH_FN(int_hash, int)
DEFINE_COMPARE_AND_EQUALS_FN(int_compare, int_equals, int)


static struct hashset set;
static hash_fn hash;
static equals_fn equals;
static int *vals;
static ssize_t size;


static void empty_setup_fixture(void **state)
{
	print_message("empty hashset\n");
	print_message("-------------\n");
}

static void teardown_fixture(void **state)
{
	print_message("\n\n");
}

static void empty_setup(void **state)
{
	static int empty_vals[] = { };

	hash = int_hash;
	equals = int_equals;
	hashset_init(&set, hash, equals, sizeof(int));
	
	vals = empty_vals;
	size = 0;
}

static void empty_teardown(void **state)
{
	hashset_deinit(&set);
}

static void test_size(void **state)
{
	assert_int_equal(hashset_size(&set), size);
}

static void test_add(void **state)
{
	int val = 31337;
	
	hashset_add(&set, &val);
	assert_int_equal(hashset_size(&set), size + 1);
	assert_true(hashset_contains(&set, &val));
	assert_int_equal(*(int *)hashset_lookup(&set, &val), val);
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_size, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_add, empty_setup, empty_teardown),
		unit_test_teardown(empty_suite, teardown_fixture),
	};
	return run_tests(tests);
}


