#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <checkutils.h>
#include <iproc/actors.h>
#include <iproc/vector.h>

iproc_actors *actors;

static void
empty_setup ()
{
    double defvector_data[] = { -1.0, -1.0, -1.0 };
    iproc_vector_view def = iproc_vector_view_array(defvector_data, 3);

    actors = iproc_actors_new(5, &def.vector);
}

static void
teardown ()
{
    iproc_actors_unref(actors);
}

START_TEST (test_empty_init)
{
    ck_assert(actors != NULL);
    ck_assert_int_eq(iproc_actors_size(actors), 5);
    ck_assert_int_eq(iproc_actors_nclass(actors), 1);
    ck_assert_int_eq(iproc_actors_dim(actors), 3);

}
END_TEST

START_TEST (test_empty_class)
{
    ck_assert_int_eq(iproc_actors_class(actors, 0), IPROC_ACTORS_DEFCLASS);
    ck_assert_int_eq(iproc_actors_class(actors, 4), IPROC_ACTORS_DEFCLASS);
}
END_TEST

START_TEST (test_empty_vector)
{
    iproc_vector *x;

    x = iproc_actors_vector(actors, 4);
    ck_assert(x != NULL);
    ck_assert_int_eq(iproc_vector_dim(x), 3);
    ck_assert_feq(iproc_vector_get(x, 0), -1.0);
    ck_assert_feq(iproc_vector_get(x, 1), -1.0);
    ck_assert_feq(iproc_vector_get(x, 2), -1.0);

    x = iproc_actors_vector(actors, 0);
    ck_assert(x != NULL);
    ck_assert_int_eq(iproc_vector_dim(x), 3);
    ck_assert_feq(iproc_vector_get(x, 0), -1.0);
    ck_assert_feq(iproc_vector_get(x, 1), -1.0);
    ck_assert_feq(iproc_vector_get(x, 2), -1.0);
}
END_TEST

START_TEST (test_empty_class_vector)
{
    iproc_vector *x = iproc_actors_class_vector(actors, IPROC_ACTORS_DEFCLASS);
    ck_assert(x != NULL);
    ck_assert_int_eq(iproc_vector_dim(x), 3);
    ck_assert_feq(iproc_vector_get(x, 0), -1.0);
    ck_assert_feq(iproc_vector_get(x, 1), -1.0);
    ck_assert_feq(iproc_vector_get(x, 2), -1.0);
}
END_TEST


static void
twoclass_setup ()
{
    double defvector_data[] = { 0.0, 2.0, 4.0 };
    iproc_vector_view x0 = iproc_vector_view_array(defvector_data + 0, 1);
    iproc_vector_view x1 = iproc_vector_view_array(defvector_data + 1, 1);
    iproc_vector_view x2 = iproc_vector_view_array(defvector_data + 2, 1);

    actors = iproc_actors_new(7, &x0.vector);
    iproc_actors_append_class(actors, &x1.vector);
    iproc_actors_append_class(actors, &x2.vector);

    iproc_actors_set(actors, 1, 2);
    iproc_actors_set(actors, 5, 1);
    iproc_actors_set(actors, 3, 2);
}

START_TEST (test_twoclass_init)
{
    ck_assert(actors != NULL);
    ck_assert_int_eq(iproc_actors_size(actors), 7);
    ck_assert_int_eq(iproc_actors_nclass(actors), 3);
    ck_assert_int_eq(iproc_actors_dim(actors), 1);

}
END_TEST

START_TEST (test_twoclass_class)
{
    ck_assert_int_eq(iproc_actors_class(actors, 1), 2);
    ck_assert_int_eq(iproc_actors_class(actors, 3), 2);
    ck_assert_int_eq(iproc_actors_class(actors, 5), 1);

    ck_assert_int_eq(iproc_actors_class(actors, 0), 0);
    ck_assert_int_eq(iproc_actors_class(actors, 2), 0);
    ck_assert_int_eq(iproc_actors_class(actors, 4), 0);
    ck_assert_int_eq(iproc_actors_class(actors, 6), 0);
}
END_TEST

START_TEST (test_twoclass_mul_notrans)
{
    double x_data[] = { 3.0 };
    double y_data[] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
    int n = 7;
    int p = 1;
    iproc_vector_view x = iproc_vector_view_array(x_data, p);
    iproc_vector_view y = iproc_vector_view_array(y_data, n);

    iproc_actors_mul(1.0, IPROC_TRANS_NOTRANS, actors, &x.vector, 0.0, &y.vector);
    ck_assert_feq(y_data[0],  0.0);
    ck_assert_feq(y_data[1], 12.0);
    ck_assert_feq(y_data[2],  0.0);
    ck_assert_feq(y_data[3], 12.0);
    ck_assert_feq(y_data[4],  0.0);
    ck_assert_feq(y_data[5],  6.0);
    ck_assert_feq(y_data[6],  0.0);

    iproc_actors_mul(2.0, IPROC_TRANS_NOTRANS, actors, &x.vector, 1.0, &y.vector);
    ck_assert_feq(y_data[0],  0.0);
    ck_assert_feq(y_data[1], 36.0);
    ck_assert_feq(y_data[2],  0.0);
    ck_assert_feq(y_data[3], 36.0);
    ck_assert_feq(y_data[4],  0.0);
    ck_assert_feq(y_data[5], 18.0);
    ck_assert_feq(y_data[6],  0.0);

    iproc_actors_mul(-1.0, IPROC_TRANS_NOTRANS, actors, &x.vector, 0.5, &y.vector);
    ck_assert_feq(y_data[0],  0.0);
    ck_assert_feq(y_data[1],  6.0);
    ck_assert_feq(y_data[2],  0.0);
    ck_assert_feq(y_data[3],  6.0);
    ck_assert_feq(y_data[4],  0.0);
    ck_assert_feq(y_data[5],  3.0);
    ck_assert_feq(y_data[6],  0.0);
}
END_TEST

START_TEST (test_twoclass_mul_trans)
{
    double x_data[] = { 1.0, 0.0, 1.0, 2.0, 1.0, -1.0, 1.0 };
    double y_data[] = { 2.0 };
    int n = 7;
    int p = 1;
    iproc_vector_view x = iproc_vector_view_array(x_data, n);
    iproc_vector_view y = iproc_vector_view_array(y_data, p);

    iproc_actors_mul(1.0, IPROC_TRANS_TRANS, actors, &x.vector, 0.0, &y.vector);
    ck_assert_feq(y_data[0], 6.0);

    iproc_actors_mul(2.0, IPROC_TRANS_TRANS, actors, &x.vector, 1.0, &y.vector);
    ck_assert_feq(y_data[0], 18.0);

    iproc_actors_mul(-1.0, IPROC_TRANS_TRANS, actors, &x.vector, 0.5, &y.vector);
    ck_assert_feq(y_data[0], 3.0);
}
END_TEST

static TCase *
empty_tcase ()
{
    TCase *tc = tcase_create("empty actors");
    tcase_add_checked_fixture(tc, empty_setup, teardown);
    tcase_add_test(tc, test_empty_init);
    tcase_add_test(tc, test_empty_class);
    tcase_add_test(tc, test_empty_vector);
    tcase_add_test(tc, test_empty_class_vector);
    return tc;
}

static TCase *
twoclass_tcase ()
{
    TCase *tc = tcase_create("actors with two classes");
    tcase_add_checked_fixture(tc, twoclass_setup, teardown);
    tcase_add_test(tc, test_twoclass_init);
    tcase_add_test(tc, test_twoclass_class);
    tcase_add_test(tc, test_twoclass_mul_notrans);
    tcase_add_test(tc, test_twoclass_mul_trans);
    return tc;
}


Suite *
actors_suite ()
{
    Suite *s = suite_create("actors");
    suite_add_tcase(s, empty_tcase());
    suite_add_tcase(s, twoclass_tcase());
    return s;
}
