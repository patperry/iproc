#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <checkutils.h>
#include <iproc/svector.h>

iproc_svector *x;

static void
zero5_setup ()
{
    x = iproc_svector_new(5);
}

static void
teardown ()
{
    iproc_svector_unref(x);
}

START_TEST (test_z5_init)
{
    ck_assert(x != NULL);
    ck_assert_int_eq(iproc_svector_dim(x), 5);
    ck_assert_int_eq(iproc_svector_nnz(x), 0);
}
END_TEST

START_TEST (test_z5_get)
{
    int64_t i;

    for (i = 0; i < 5; i++) {
        ck_assert_feq(iproc_svector_get(x, i), 0.0);
    }
}
END_TEST

START_TEST (test_z5_set)
{
    iproc_svector_set(x, 4, 1.23);

    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), 0.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);

    iproc_svector_set(x, 2, -222.0);

    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), -222.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);

    iproc_svector_set(x, 3, 3.14);
    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), -222.0);
    ck_assert_feq(iproc_svector_get(x, 3), 3.14);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);

}
END_TEST

START_TEST (test_z5_inc)
{
    iproc_svector_inc(x, 4, 1.23);

    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), 0.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);

    iproc_svector_inc(x, 2, -222.0);
    iproc_svector_inc(x, 2, -222.0);

    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), -444.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);

    iproc_svector_inc(x, 3, 3.14);
    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), -444.0);
    ck_assert_feq(iproc_svector_get(x, 3), 3.14);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);

    iproc_svector_inc(x, 2, -222.0);
    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 0.0);
    ck_assert_feq(iproc_svector_get(x, 2), -666.0);
    ck_assert_feq(iproc_svector_get(x, 3), 3.14);
    ck_assert_feq(iproc_svector_get(x, 4), 1.23);
}
END_TEST

static void
basis4_setup ()
{
    x = iproc_svector_new(4);
    iproc_svector_set(x, 1, 1.0);
}

START_TEST (test_b4_init)
{
    ck_assert(x != NULL);
    ck_assert_int_eq(iproc_svector_dim(x), 4);
    ck_assert_int_eq(iproc_svector_nnz(x), 1);
}
END_TEST

START_TEST (test_b4_get)
{
    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 1.0);
    ck_assert_feq(iproc_svector_get(x, 2), 0.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);
}
END_TEST

START_TEST (test_b4_set)
{
    iproc_svector_set(x, 1, 1.23);
    ck_assert_feq(iproc_svector_get(x, 0), 0.0);
    ck_assert_feq(iproc_svector_get(x, 1), 1.23);
    ck_assert_feq(iproc_svector_get(x, 2), 0.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);


    iproc_svector_set(x, 0, 23);
    ck_assert_feq(iproc_svector_get(x, 0), 23);
    ck_assert_feq(iproc_svector_get(x, 1), 1.23);
    ck_assert_feq(iproc_svector_get(x, 2), 0.0);
    ck_assert_feq(iproc_svector_get(x, 3), 0.0);
}
END_TEST

START_TEST (test_b4_sdot)
{
    double elems[] = { 1.0, 2.0, 3.0, 4.0 };
    iproc_vector_view y = iproc_vector_view_array(elems, 4);

    double dot = iproc_vector_sdot(&y.vector, x);
    ck_assert_feq(dot, 2.0);
}
END_TEST


START_TEST (test_b4_sacc)
{
    double elems[] = { 1.0, 2.0, 3.0, 4.0 };
    iproc_vector_view y = iproc_vector_view_array(elems, 4);

    iproc_vector_sacc(&y.vector, -4.2, x);
    ck_assert_feq(elems[0],  1.0);
    ck_assert_feq(elems[1], -2.2);
    ck_assert_feq(elems[2],  3.0);
    ck_assert_feq(elems[3],  4.0);
}
END_TEST


static TCase *
zero5_tcase ()
{
    TCase *tc = tcase_create("zero 5-svector");
    tcase_add_checked_fixture(tc, zero5_setup, teardown);
    tcase_add_test(tc, test_z5_init);
    tcase_add_test(tc, test_z5_get);
    tcase_add_test(tc, test_z5_set);
    tcase_add_test(tc, test_z5_inc);
    return tc;
}

static TCase *
basis4_tcase ()
{
    TCase *tc = tcase_create("basis 4-svector");
    tcase_add_checked_fixture(tc, basis4_setup, teardown);
    tcase_add_test(tc, test_b4_init);
    tcase_add_test(tc, test_b4_get);
    tcase_add_test(tc, test_b4_set);
    tcase_add_test(tc, test_b4_sdot);
    tcase_add_test(tc, test_b4_sacc);
    return tc;
}


Suite *
svector_suite ()
{
    Suite *s = suite_create("svector");
    suite_add_tcase(s, zero5_tcase());
    suite_add_tcase(s, basis4_tcase());
    return s;
}
