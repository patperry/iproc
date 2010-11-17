
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include <libla/vector.h>


LAVector *v;
la_size n;
double *elems_v;

START_TEST (test_new)
{
    fail_unless(v != NULL,
                "Vector is non-null");
    fail_unless(la_vector_dim(v) == n,
                "Dimension not set correctly on creation");
}
END_TEST

START_TEST (test_get)
{
    int i;
    for (i = 0; i < n; i++) {
        fail_unless(la_vector_get(v, i) == elems_v[i],
                    "Element value %d not retreived correctly", i);
    }
}
END_TEST

START_TEST (test_set_all)
{
    int i;
    la_vector_set_all(v, 1234);
    
    for (i = 0; i < n; i++) {
        fail_unless(la_vector_get(v, i) == 1234,
                    "Element %d not initialized", i);
    }
}
END_TEST

START_TEST (test_set_basis)
{
    int i;
    
    la_vector_set_basis(v, _i);

    fail_unless(la_vector_get(v, _i) == 1);
    for (i = 0; i < n; i++) {
        if (i == _i)
            continue;

        fail_unless(la_vector_get(v, i) == 0);
    }
}
END_TEST

START_TEST (test_subvector)
{
    int i;
    LAVectorView v3 = la_vector_subvector(v, 1, 3);
    la_vector_set_all(&v3.vector, 100);
    fail_unless(la_vector_get(v, 0) == elems_v[0]);
    for (i = 1; i < 4; i++) {
        fail_unless(la_vector_get(v, i) == 100);
    }
    fail_unless(la_vector_get(v, 4) == elems_v[4]);
}
END_TEST


void
v0_setup ()
{
    n = 0;
    v = la_vector_new(0);
}

void
v0_teardown()
{
    la_vector_free(v);
}

TCase *
v0_tcase ()
{
    TCase *tc = tcase_create("vector with dimension 0");
    tcase_add_checked_fixture(tc, v0_setup, v0_teardown);
    tcase_add_test(tc, test_new);
    tcase_add_test(tc, test_get);
    tcase_add_test(tc, test_set_all);
    return tc;
}


void
v5_setup ()
{
    double elems[5] = { 1.12, 0, -23, 0, -0.0001 };
    int i;

    n = 5;
    v = la_vector_new(5);
    elems_v = malloc(n * sizeof(double));
    memcpy(elems_v, elems, n * sizeof(double));
    memcpy(la_vector_ptr(v, 0), elems, n * sizeof(double));
}

void
v5_teardown ()
{
    free(elems_v);
    la_vector_free(v);
}

TCase *
v5_tcase ()
{
    TCase *tc = tcase_create("vector with dimension 5");
    tcase_add_checked_fixture(tc, v5_setup, v5_teardown);
    tcase_add_test(tc, test_new);
    tcase_add_test(tc, test_get);
    tcase_add_test(tc, test_set_all);
    tcase_add_loop_test(tc, test_set_basis, 0, 5);
    tcase_add_test(tc, test_subvector);
    return tc;
}

Suite *
vector_suite ()
{
    Suite *s = suite_create("vector");
    suite_add_tcase(s, v0_tcase());
    suite_add_tcase(s, v5_tcase());
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
