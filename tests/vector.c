
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <check.h>
#include <checkutils.h>
#include <libla/vector.h>
#include <libla/ieee754.h>

LAVector *v;
la_size n;
double *elems_v;

LAVector *v1, *v2;
double *elems_v1;
double *elems_v2;

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
        ck_assert_feq(la_vector_get(v, i), elems_v[i]);
    }
}
END_TEST

START_TEST (test_set_all)
{
    int i;
    la_vector_set_all(v, 1234);

    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v, i), 1234);
    }
}
END_TEST

START_TEST (test_reverse)
{
    int i;

    la_vector_reverse(v);
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v, i), elems_v[n - i - 1]);
    }
}
END_TEST

START_TEST (test_scale)
{
    int i;
    la_vector_scale(v, 2.0);
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v, i), 2.0 * elems_v[i]);
    }
}
END_TEST

START_TEST (test_shift)
{
    int i;
    la_vector_shift(v, 1.2);
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v, i), 1.2 + elems_v[i]);
    }
}
END_TEST


TCase *
vector_tcase (char *desc, void (*setup)(), void (*teardown)())
{
    TCase *tc = tcase_create(desc);
    tcase_add_checked_fixture(tc, setup, teardown);
    tcase_add_test(tc, test_new);
    tcase_add_test(tc, test_get);
    tcase_add_test(tc, test_set_all);
    tcase_add_test(tc, test_reverse);
    tcase_add_test(tc, test_scale);
    tcase_add_test(tc, test_shift);

    return tc;
}


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
    TCase *tc = vector_tcase("vector with dimension 0",
                             v0_setup,
                             v0_teardown);
    return tc;
}


void
v5_setup ()
{
    double elems[5] = { 1.23, 0, -21, 0, -0.0001 };

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

START_TEST (v5_test_set_basis)
{
    int i;

    la_vector_set_basis(v, 2);

    ck_assert_feq(la_vector_get(v, 2), 1);
    for (i = 0; i < n; i++) {
        if (i == 2)
            continue;

        ck_assert_feq(la_vector_get(v, i), 0);
    }
}
END_TEST

START_TEST (v5_test_subvector)
{
    int i;
    LAVectorView v3 = la_vector_subvector(v, 1, 3);
    la_vector_set_all(&v3.vector, 100);
    ck_assert_feq(la_vector_get(v, 0), elems_v[0]);
    for (i = 1; i < 4; i++) {
        ck_assert_feq(la_vector_get(v, i), 100);
    }
    ck_assert_feq(la_vector_get(v, 4), elems_v[4]);
}
END_TEST

START_TEST (v5_test_swap_elems)
{
    la_vector_swap_elems(v, 0, 4);
    ck_assert_feq(la_vector_get(v, 0), elems_v[4]);
    ck_assert_feq(la_vector_get(v, 4), elems_v[0]);
}
END_TEST

TCase *
v5_tcase ()
{
    TCase *tc = vector_tcase("vector with dimension 5",
                             v5_setup,
                             v5_teardown);

    tcase_add_test(tc, v5_test_set_basis);
    tcase_add_test(tc, v5_test_subvector);
    tcase_add_test(tc, v5_test_swap_elems);
    return tc;
}

START_TEST (test_add)
{
    int i;
    la_vector_add(v1, v2);

    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v1, i), elems_v1[i] + elems_v2[i]);
    }
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v2, i), elems_v2[i]);
    }
}
END_TEST


START_TEST (test_sub)
{
    int i;
    la_vector_sub(v1, v2);

    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v1, i), elems_v1[i] - elems_v2[i]);
    }
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v2, i), elems_v2[i]);
    }
}
END_TEST


START_TEST (test_mul)
{
    int i;
    la_vector_mul(v1, v2);

    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v1, i), elems_v1[i] * elems_v2[i]);
    }
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v2, i), elems_v2[i]);
    }
}
END_TEST


START_TEST (test_div)
{
    int i;
    la_vector_div(v1, v2);

    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v1, i), elems_v1[i] / elems_v2[i]);
    }
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v2, i), elems_v2[i]);
    }
}
END_TEST

START_TEST (test_acc)
{
    int i;
    la_vector_acc(v1, 3.14, v2);

    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v1, i), elems_v1[i] + 3.14 * elems_v2[i]);
    }
    for (i = 0; i < n; i++) {
        ck_assert_feq(la_vector_get(v2, i), elems_v2[i]);
    }
}
END_TEST



TCase *
vector_pair_tcase (char *desc, void (*setup)(), void (*teardown)())
{
    TCase *tc = tcase_create(desc);
    tcase_add_checked_fixture(tc, setup, teardown);
    tcase_add_test(tc, test_add);
    tcase_add_test(tc, test_sub);
    tcase_add_test(tc, test_mul);
    tcase_add_test(tc, test_div);
    tcase_add_test(tc, test_acc);
    return tc;
}

void
v0_pair_setup ()
{
    n = 0;
    v1 = la_vector_new(n);
    v2 = la_vector_new(n);
}

void
v0_pair_teardown ()
{
    la_vector_free(v1);
    la_vector_free(v2);
}

TCase *
v0_pair_tcase ()
{
    TCase *tc = vector_pair_tcase("two vectors with dimension 0",
                                  v0_pair_setup,
                                  v0_pair_teardown);
    return tc;
}

void
v5_pair_setup ()
{
    double elems1[5] = { 0, -0.7, 21.0,  0.370, -1.0 };
    double elems2[5] = { 4,  0.0,  0.8,  0.007, -0.2 };

    n = 5;
    v1 = la_vector_new(n);
    v2 = la_vector_new(n);
    elems_v1 = malloc(n * sizeof(double));
    elems_v2 = malloc(n * sizeof(double));
    memcpy(elems_v1, elems1, n * sizeof(double));
    memcpy(elems_v2, elems2, n * sizeof(double));
    memcpy(la_vector_ptr(v1, 0), elems1, n * sizeof(double));
    memcpy(la_vector_ptr(v2, 0), elems2, n * sizeof(double));
}

void
v5_pair_teardown ()
{
    free(elems_v1);
    free(elems_v2);
    la_vector_free(v1);
    la_vector_free(v2);
}

TCase *
v5_pair_tcase ()
{
    TCase *tc = vector_pair_tcase("two vectors with dimension 5",
                                  v5_pair_setup,
                                  v5_pair_teardown);
    return tc;
}

Suite *
vector_suite ()
{
    Suite *s = suite_create("vector");
    suite_add_tcase(s, v0_tcase());
    suite_add_tcase(s, v5_tcase());
    suite_add_tcase(s, v0_pair_tcase());
    suite_add_tcase(s, v5_pair_tcase());
    return s;
}
