#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <cmockery.h>

#include "pqueue.h"


static iproc_pqueue *pqueue = NULL;
static ssize_t size;
static ssize_t * elts; // elements sorted in decending order

static void
empty_setup_fixture (void **state)
{
    print_message("Empty pqueue\n");
    print_message("------------\n");
}

static void
empty_setup (void **state)
{
    static ssize_t empty_elts[] = { };
    
    pqueue = iproc_pqueue_new(sizeof(ssize_t), iproc_ssize_compare);
    size = 0;
    elts = empty_elts;
}

static void
singleton_setup_fixture (void **state)
{
    print_message("Singleton pqueue\n");
    print_message("----------------\n");    
}

static void
singleton_setup (void **state)
{
    static ssize_t singleton_elts[] = { 1234 };

    pqueue = iproc_pqueue_new(sizeof(ssize_t), iproc_ssize_compare);
    elts = singleton_elts;
    size = 1;
    iproc_pqueue_push_array(pqueue, elts, size);
}

static void
sorted5_setup_fixture (void **state)
{
    print_message("Pqueue with 5 sorted elements\n");
    print_message("-----------------------------\n");    
}

static void
sorted5_setup (void **state)
{
    static ssize_t sorted5_elts[] = { 5, 4, 3, 2, 1 };
    
    pqueue = iproc_pqueue_new(sizeof(ssize_t), iproc_ssize_compare);
    elts = sorted5_elts;
    size = 5;
    iproc_pqueue_push_array(pqueue, elts, size);

}

static void
teardown (void **state)
{
    iproc_pqueue_unref(pqueue);
}

static void
teardown_fixture (void **state)
{
    print_message("\n\n");
}



static void
test_size (void **state)
{
    assert_int_equal(iproc_pqueue_size(pqueue), size);
}

static void
test_push_min_minus_one (void **state)
{
    ssize_t min = elts[size - 1];
    ssize_t min_minus_one = min - 1;
    iproc_pqueue_push(pqueue, &min_minus_one);
    assert_int_equal(iproc_pqueue_size(pqueue), size + 1);
    assert_int_equal(*(ssize_t *)iproc_pqueue_top(pqueue), elts[0]);
}

static void
test_push_min (void **state)
{
    ssize_t min = elts[size - 1];
    ssize_t elt = min;
    iproc_pqueue_push(pqueue, &elt);
    assert_int_equal(iproc_pqueue_size(pqueue), size + 1);
    assert_int_equal(*(ssize_t *)iproc_pqueue_top(pqueue), elts[0]);
}

static void
test_push_max_minus_one (void **state)
{
    ssize_t max = elts[0];
    ssize_t max_minus_one = max - 1;
    iproc_pqueue_push(pqueue, &max_minus_one);
    assert_int_equal(iproc_pqueue_size(pqueue), size + 1);
    assert_int_equal(*(ssize_t *)iproc_pqueue_top(pqueue), max);
}

static void
test_push_max (void **state)
{
    ssize_t max = elts[0];
    ssize_t elt = max;
    iproc_pqueue_push(pqueue, &elt);
    assert_int_equal(iproc_pqueue_size(pqueue), size + 1);
    assert_int_equal(*(ssize_t *)iproc_pqueue_top(pqueue), max);
}

static void
test_push_max_plus_one (void **state)
{
    ssize_t max = (size ? elts[0] : 0);
    ssize_t max_plus_one = max + 1;
    iproc_pqueue_push(pqueue, &max_plus_one);
    assert_int_equal(iproc_pqueue_size(pqueue), size + 1);
    assert_int_equal(*(ssize_t *)iproc_pqueue_top(pqueue), max + 1);
}

int
main(int argc, char **argv)
{
    UnitTest tests[] = {
        unit_test_setup(empty_suite, empty_setup_fixture),
        unit_test_setup_teardown(test_size, empty_setup, teardown),
        unit_test_teardown(empty_suite, teardown_fixture),

        unit_test_setup(singleton_suite, singleton_setup_fixture),
        unit_test_setup_teardown(test_size,               singleton_setup, teardown),
        unit_test_setup_teardown(test_push_min,           singleton_setup, teardown),
        unit_test_setup_teardown(test_push_min_minus_one, singleton_setup, teardown),
        unit_test_setup_teardown(test_push_max,           singleton_setup, teardown),
        unit_test_setup_teardown(test_push_max_minus_one, singleton_setup, teardown),
        unit_test_setup_teardown(test_push_max_plus_one,  singleton_setup, teardown),
        unit_test_teardown(singleton_suite, teardown_fixture),
        
        unit_test_setup(sorted5_suite, sorted5_setup_fixture),
        unit_test_setup_teardown(test_size,               sorted5_setup, teardown),
        unit_test_setup_teardown(test_push_min,           sorted5_setup, teardown),
        unit_test_setup_teardown(test_push_min_minus_one, sorted5_setup, teardown),
        unit_test_setup_teardown(test_push_max,           sorted5_setup, teardown),
        unit_test_setup_teardown(test_push_max_minus_one, sorted5_setup, teardown),
        unit_test_setup_teardown(test_push_max_plus_one,  sorted5_setup, teardown),
        unit_test_teardown(sorted5_suite, teardown_fixture)
    };
    return run_tests(tests);
}
