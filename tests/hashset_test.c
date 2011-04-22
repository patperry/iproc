#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "hashset.h"

DEFINE_HASH_FN(int_hash, int)
DEFINE_COMPARE_AND_EQUALS_FN(int_compare, int_equals, int)

static uint32_t int_bad_hash(const void *x)
{
	return 1337;
}

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
	static int *empty_vals = NULL;

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

static void big_setup_fixture(void **state)
{
	print_message("big hashset\n");
	print_message("-----------\n");
}

static void big_setup(void **state)
{
	empty_setup(state);
	
	size = 555;
	vals = malloc(size * sizeof(*vals));
	ssize_t i;
	
	for (i = 0; i < size; i++) {
		vals[i] = (int)i;
	}

	hashset_add_all(&set, vals, size);
}

static void big_teardown(void **state)
{
	free(vals);
	empty_teardown(state);
}

static void big_bad_setup_fixture(void **state)
{
	print_message("big hashset (bad hash)\n");
	print_message("----------------------\n");
}

static void big_bad_setup(void **state)
{
	hash = int_bad_hash;
	equals = int_equals;
	hashset_init(&set, hash, equals, sizeof(int));
	
	size = 151;
	vals = malloc(size * sizeof(*vals));
	ssize_t i;
	
	for (i = 0; i < size; i++) {
		vals[i] = (int)i;
	}
	
	hashset_add_all(&set, vals, size);
}

static void big_bad_teardown(void **state)
{
	free(vals);
	empty_teardown(state);
}


static void test_size(void **state)
{
	assert_int_equal(hashset_size(&set), size);
	if (size == 0) {
		assert_true(hashset_empty(&set));
	} else {
		assert_false(hashset_empty(&set));
	}
}

static void test_clear(void **state)
{
	hashset_clear(&set);
	assert_true(hashset_empty(&set));
}

static void test_lookup(void **state)
{
	ssize_t i;
	const int *val;
	
	for (i = 0; i < size; i++) {
		assert_true(hashset_contains(&set, &vals[i]));
		val = hashset_lookup(&set, &vals[i]);
		assert_true(val);
		assert_int_equal(*val, vals[i]);
	}
}

static void test_add(void **state)
{
	int val = 31337;
	
	hashset_add(&set, &val);
	assert_int_equal(hashset_size(&set), size + 1);
	assert_true(hashset_contains(&set, &val));
	assert_int_equal(*(int *)hashset_lookup(&set, &val), val);
}

static void test_add_existing(void **state)
{
	int val = 88888;
	
	hashset_add(&set, &val);
	hashset_add(&set, &val);
	assert_int_equal(hashset_size(&set), size + 1);
	assert_true(hashset_contains(&set, &val));
	assert_int_equal(*(int *)hashset_lookup(&set, &val), val);
}

static void test_remove(void **state)
{
	int val = -1;
	
	hashset_add(&set, &val);
	hashset_remove(&set, &val);
	assert_int_equal(hashset_size(&set), size);	
	assert_false(hashset_contains(&set, &val));
	assert_false(hashset_lookup(&set, &val));
}

static void test_remove_hard(void **state)
{
	ssize_t i, j;
	
	for (i = 0; i < size; i++) {
		hashset_remove(&set, &vals[i]);
		assert_int_equal(hashset_size(&set), size - i - 1);
		for (j = 0; j <= i; j++) {
			assert_false(hashset_contains(&set, &vals[j]));
		}
		for (; j < size; j++) {
			assert_true(hashset_contains(&set, &vals[j]));
			assert_int_equal(*(int *)hashset_lookup(&set, &vals[j]), vals[j]);
		}
	}
	assert_true(hashset_empty(&set));
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_size, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_clear, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_lookup, empty_setup, empty_teardown),		
		unit_test_setup_teardown(test_add, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_add_existing, empty_setup, empty_teardown),
		unit_test_setup_teardown(test_remove, empty_setup, empty_teardown),		
		unit_test_teardown(empty_suite, teardown_fixture),

		unit_test_setup(big_suite, big_setup_fixture),
		unit_test_setup_teardown(test_size, big_setup, big_teardown),
		unit_test_setup_teardown(test_clear, big_setup, big_teardown),		
		unit_test_setup_teardown(test_lookup, big_setup, big_teardown),		
		unit_test_setup_teardown(test_add, big_setup, big_teardown),
		unit_test_setup_teardown(test_add_existing, big_setup, big_teardown),
		unit_test_setup_teardown(test_remove, big_setup, big_teardown),		
		unit_test_setup_teardown(test_remove_hard, big_setup, big_teardown),
		unit_test_teardown(big_suite, teardown_fixture),

		unit_test_setup(big_bad_suite, big_bad_setup_fixture),
		unit_test_setup_teardown(test_size, big_bad_setup, big_bad_teardown),
		unit_test_setup_teardown(test_clear, big_bad_setup, big_bad_teardown),				
		unit_test_setup_teardown(test_lookup, big_bad_setup, big_bad_teardown),		
		unit_test_setup_teardown(test_add, big_bad_setup, big_teardown),
		unit_test_setup_teardown(test_add_existing, big_bad_setup, big_bad_teardown),
		unit_test_setup_teardown(test_remove, big_bad_setup, big_bad_teardown),		
		unit_test_setup_teardown(test_remove_hard, big_bad_setup, big_bad_teardown),
		unit_test_teardown(big_bad_suite, teardown_fixture),
	};
	return run_tests(tests);
}


