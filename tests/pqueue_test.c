#include "port.h"

#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "pqueue.h"

DEFINE_COMPARE_AND_EQUALS_FN(ssize_compare, ssize_equals, ssize_t)

static struct pqueue pqueue;
static ssize_t size;
static ssize_t *elts;		// elements sorted in decending order

static void empty_setup_fixture(void **state)
{
	print_message("Empty pqueue\n");
	print_message("------------\n");
}

static void empty_setup(void **state)
{
	static ssize_t *empty_elts = NULL;

	pqueue_init(&pqueue, ssize_compare, NULL, sizeof(ssize_t));
	size = 0;
	elts = empty_elts;
}

static void singleton_setup_fixture(void **state)
{
	print_message("Singleton pqueue\n");
	print_message("----------------\n");
}

static void singleton_setup(void **state)
{
	static ssize_t singleton_elts[] = { 1234 };

	pqueue_init(&pqueue, ssize_compare, NULL, sizeof(ssize_t));
	elts = singleton_elts;
	size = 1;
	pqueue_push_all(&pqueue, elts, size);
}

static void sorted5_setup_fixture(void **state)
{
	print_message("Pqueue with 5 sorted elements\n");
	print_message("-----------------------------\n");
}

static void sorted5_setup(void **state)
{
	static ssize_t sorted5_elts[] = { 5, 4, 3, 2, 1 };

	pqueue_init(&pqueue, ssize_compare, NULL, sizeof(ssize_t));
	elts = sorted5_elts;
	size = 5;
	pqueue_push_all(&pqueue, elts, size);

}

static void unsorted7_setup_fixture(void **state)
{
	print_message("Pqueue with 7 unsorted elements\n");
	print_message("-------------------------------\n");
}

static void unsorted7_setup(void **state)
{
	static ssize_t sorted7_elts[] = { 7, 6, 5, 4, 3, 2, 1 };
	ssize_t unsorted7_elts[] = { 2, 1, 3, 4, 7, 6, 5 };

	pqueue_init(&pqueue, ssize_compare, NULL, sizeof(ssize_t));
	elts = sorted7_elts;
	size = 7;
	pqueue_push_all(&pqueue, unsorted7_elts, size);

}

static void teardown(void **state)
{
	pqueue_deinit(&pqueue);
}

static void teardown_fixture(void **state)
{
	print_message("\n\n");
}

static void test_size(void **state)
{
	assert_int_equal(pqueue_size(&pqueue), size);
}

static void test_push_min_minus_one(void **state)
{
	ssize_t min = elts[size - 1];
	ssize_t min_minus_one = min - 1;
	pqueue_push(&pqueue, &min_minus_one);
	assert_int_equal(pqueue_size(&pqueue), size + 1);
	assert_int_equal(*(ssize_t *)pqueue_top(&pqueue), elts[0]);
}

static void test_push_min(void **state)
{
	ssize_t min = elts[size - 1];
	ssize_t elt = min;
	pqueue_push(&pqueue, &elt);
	assert_int_equal(pqueue_size(&pqueue), size + 1);
	assert_int_equal(*(ssize_t *)pqueue_top(&pqueue), elts[0]);
}

static void test_push_max_minus_one(void **state)
{
	ssize_t max = elts[0];
	ssize_t max_minus_one = max - 1;
	pqueue_push(&pqueue, &max_minus_one);
	assert_int_equal(pqueue_size(&pqueue), size + 1);
	assert_int_equal(*(ssize_t *)pqueue_top(&pqueue), max);
}

static void test_push_max(void **state)
{
	ssize_t max = elts[0];
	ssize_t elt = max;
	pqueue_push(&pqueue, &elt);
	assert_int_equal(pqueue_size(&pqueue), size + 1);
	assert_int_equal(*(ssize_t *)pqueue_top(&pqueue), max);
}

static void test_push_max_plus_one(void **state)
{
	ssize_t max = (size ? elts[0] : 0);
	ssize_t max_plus_one = max + 1;
	pqueue_push(&pqueue, &max_plus_one);
	assert_int_equal(pqueue_size(&pqueue), size + 1);
	assert_int_equal(*(ssize_t *)pqueue_top(&pqueue), max + 1);
}

static void test_push_existing(void **state)
{
	ssize_t i, j;
	ssize_t top;
	struct pqueue pq;
	ssize_t elt;

	for (i = 0; i < size; i++) {
		pqueue_init_copy(&pq, &pqueue);
		elt = elts[i];

		pqueue_push(&pq, &elt);

		assert_int_equal(pqueue_size(&pq), size + 1);

		for (j = 0; j < size + 1; j++) {
			top = *(ssize_t *)pqueue_top(&pq);
			pqueue_pop(&pq);
			if (j <= i) {
				assert_int_equal(top, elts[j]);
			} else {
				assert_int_equal(top, elts[j - 1]);
			}
		}

		pqueue_deinit(&pq);
	}
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_size, empty_setup, teardown),
		unit_test_teardown(empty_suite, teardown_fixture),

		unit_test_setup(singleton_suite, singleton_setup_fixture),
		unit_test_setup_teardown(test_size, singleton_setup, teardown),
		unit_test_setup_teardown(test_push_min, singleton_setup,
					 teardown),
		unit_test_setup_teardown(test_push_min_minus_one,
					 singleton_setup, teardown),
		unit_test_setup_teardown(test_push_max, singleton_setup,
					 teardown),
		unit_test_setup_teardown(test_push_max_minus_one,
					 singleton_setup, teardown),
		unit_test_setup_teardown(test_push_max_plus_one,
					 singleton_setup, teardown),
		unit_test_setup_teardown(test_push_existing, singleton_setup,
					 teardown),
		unit_test_teardown(singleton_suite, teardown_fixture),

		unit_test_setup(sorted5_suite, sorted5_setup_fixture),
		unit_test_setup_teardown(test_size, sorted5_setup, teardown),
		unit_test_setup_teardown(test_push_min, sorted5_setup,
					 teardown),
		unit_test_setup_teardown(test_push_min_minus_one, sorted5_setup,
					 teardown),
		unit_test_setup_teardown(test_push_max, sorted5_setup,
					 teardown),
		unit_test_setup_teardown(test_push_max_minus_one, sorted5_setup,
					 teardown),
		unit_test_setup_teardown(test_push_max_plus_one, sorted5_setup,
					 teardown),
		unit_test_setup_teardown(test_push_existing, sorted5_setup,
					 teardown),
		unit_test_teardown(sorted5_suite, teardown_fixture),

		unit_test_setup(unsorted7_suite, unsorted7_setup_fixture),
		unit_test_setup_teardown(test_size, unsorted7_setup, teardown),
		unit_test_setup_teardown(test_push_min, unsorted7_setup,
					 teardown),
		unit_test_setup_teardown(test_push_min_minus_one,
					 unsorted7_setup, teardown),
		unit_test_setup_teardown(test_push_max, unsorted7_setup,
					 teardown),
		unit_test_setup_teardown(test_push_max_minus_one,
					 unsorted7_setup, teardown),
		unit_test_setup_teardown(test_push_max_plus_one,
					 unsorted7_setup, teardown),
		unit_test_setup_teardown(test_push_existing, unsorted7_setup,
					 teardown),
		unit_test_teardown(unsorted7_suite, teardown_fixture),

	};
	return run_tests(tests);
}
