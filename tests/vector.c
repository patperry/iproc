
#include <stdlib.h>
#include <check.h>
#include <libla/vector.h>

LAVector *v5;
const double v5_elems[5] = { 1.12, 0, -23, 0, -0.0001 };

void
v5_setup ()
{
    int i;
    v5 = la_vector_new(5);

    for (i = 0; i < 5; i++)
        la_vector_set(v5, i, v5_elems[i]);
}

void
v5_teardown ()
{
    la_vector_free(v5);
}

START_TEST (test_vector_new)
{
    fail_unless(la_vector_dim(v5) == 5,
                "Dimension not set correctly on creation");
}
END_TEST

START_TEST (test_vector_get_set)
{
    fail_unless(la_vector_get(v5, _i) == v5_elems[_i],
                "Element value not set or retreived correctly");
}
END_TEST

Suite *
vector_suite ()
{
    Suite *s = suite_create("vector");

    /* Core test case */
    TCase *tc_core = tcase_create("Core");
    tcase_add_checked_fixture(tc_core, v5_setup, v5_teardown);
    tcase_add_test(tc_core, test_vector_new);
    tcase_add_loop_test(tc_core, test_vector_get_set, 0, 5);
    suite_add_tcase(s, tc_core);

    return s;
}


int
main ()
{
    int nfail;
    Suite *s = vector_suite();
    SRunner *sr = srunner_create(s);
    srunner_run_all(sr, CK_ENV);
    nfail = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (nfail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
