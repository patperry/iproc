#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/array.h>

iproc_array *array;


static void
int_setup ()
{
    array = iproc_array_new(sizeof(int));
}

static void
charz_setup ()
{
    array = iproc_array_new(sizeof(char));
    iproc_array_append(array, "z");
}

static void
int246_setup ()
{
    int elems[] = { 2, 4, 6 };
    array = iproc_array_new(sizeof(int));
    iproc_array_append(array, elems + 0);
    iproc_array_append(array, elems + 1);
    iproc_array_append(array, elems + 2);
}

static void
teardown ()
{
    iproc_array_free(array);
}

START_TEST (test_int_init)
{
    ck_assert(array != NULL);
    ck_assert_int_eq(iproc_array_size(array), 0);
    ck_assert_int_eq(iproc_array_elem_size(array), sizeof(int));
}
END_TEST

START_TEST (test_charz_init)
{
    ck_assert(array != NULL);
    ck_assert_int_eq(iproc_array_size(array), 1);
    ck_assert_int_eq(iproc_array_elem_size(array), sizeof(char));
    ck_assert_int_eq(iproc_array_index(array, char, 0), 'z');
}
END_TEST

START_TEST (test_int246_init)
{
    ck_assert(array != NULL);
    ck_assert_int_eq(iproc_array_size(array), 3);
    ck_assert_int_eq(iproc_array_elem_size(array), sizeof(int));
    ck_assert_int_eq(iproc_array_index(array, int, 0), 2);
    ck_assert_int_eq(iproc_array_index(array, int, 1), 4);
    ck_assert_int_eq(iproc_array_index(array, int, 2), 6);
}
END_TEST

START_TEST (test_int_append)
{
    int e = 1;
    iproc_array_append(array, &e);
    ck_assert_int_eq(iproc_array_size(array), 1);
    ck_assert_int_eq(iproc_array_index(array, int, 0), 1);
}
END_TEST

START_TEST (test_int_prepend)
{
    int e = 8;
    iproc_array_prepend(array, &e);
    ck_assert_int_eq(iproc_array_size(array), 1);
    ck_assert_int_eq(iproc_array_index(array, int, 0), 8);
}
END_TEST

START_TEST (test_int_insert)
{
    int e = 18;
    iproc_array_insert(array, 1, &e);
    ck_assert_int_eq(iproc_array_size(array), 2);
    ck_assert_int_eq(iproc_array_index(array, int, 0), 0);
    ck_assert_int_eq(iproc_array_index(array, int, 1), 18);
}
END_TEST

START_TEST (test_charz_append)
{
    iproc_array_append(array, "y");
    ck_assert_int_eq(iproc_array_size(array), 2);
    ck_assert_int_eq(iproc_array_index(array, char, 0), 'z');
    ck_assert_int_eq(iproc_array_index(array, char, 1), 'y');
}
END_TEST

START_TEST (test_charz_prepend)
{
    iproc_array_prepend(array, "c");
    ck_assert_int_eq(iproc_array_size(array), 2);
    ck_assert_int_eq(iproc_array_index(array, char, 0), 'c');
    ck_assert_int_eq(iproc_array_index(array, char, 1), 'z');
}
END_TEST

START_TEST (test_charz_insert)
{
    iproc_array_insert(array, 3, "x");
    ck_assert_int_eq(iproc_array_size(array), 4);
    ck_assert_int_eq(iproc_array_index(array, char, 0), 'z');
    ck_assert_int_eq(iproc_array_index(array, char, 1), '\0');
    ck_assert_int_eq(iproc_array_index(array, char, 2), '\0');
    ck_assert_int_eq(iproc_array_index(array, char, 3), 'x');
}
END_TEST

START_TEST (test_charz_remove)
{
    iproc_array_remove(array, 0);
    ck_assert_int_eq(iproc_array_size(array), 0);
}
END_TEST

START_TEST (test_int246_append)
{
    int e = -7;
    iproc_array_append(array, &e);
    ck_assert_int_eq(iproc_array_size(array), 4);
    ck_assert_int_eq(iproc_array_index(array, int, 0),  2);
    ck_assert_int_eq(iproc_array_index(array, int, 1),  4);
    ck_assert_int_eq(iproc_array_index(array, int, 2),  6);
    ck_assert_int_eq(iproc_array_index(array, int, 3), -7);
}
END_TEST

START_TEST (test_int246_prepend)
{
    int e = 100;
    iproc_array_prepend(array, &e);
    ck_assert_int_eq(iproc_array_size(array), 4);
    ck_assert_int_eq(iproc_array_index(array, int, 0),  100);
    ck_assert_int_eq(iproc_array_index(array, int, 1),  2);
    ck_assert_int_eq(iproc_array_index(array, int, 2),  4);
    ck_assert_int_eq(iproc_array_index(array, int, 3),  6);
}
END_TEST

START_TEST (test_int246_insert)
{
    int e = 123;
    iproc_array_insert(array, 1, &e);
    ck_assert_int_eq(iproc_array_size(array), 4);
    ck_assert_int_eq(iproc_array_index(array, int, 0),   2);
    ck_assert_int_eq(iproc_array_index(array, int, 1), 123);
    ck_assert_int_eq(iproc_array_index(array, int, 2),   4);
    ck_assert_int_eq(iproc_array_index(array, int, 3),   6);
}
END_TEST

START_TEST (test_int246_remove)
{
    iproc_array_remove(array, 1);
    ck_assert_int_eq(iproc_array_size(array), 2);
    ck_assert_int_eq(iproc_array_index(array, int, 0), 2);
    ck_assert_int_eq(iproc_array_index(array, int, 1), 6);
}
END_TEST


static TCase *
int_tcase ()
{
    TCase *tc = tcase_create("int { }");
    tcase_add_checked_fixture(tc, int_setup, teardown);
    tcase_add_test(tc, test_int_init);
    tcase_add_test(tc, test_int_append);
    tcase_add_test(tc, test_int_prepend);
    tcase_add_test(tc, test_int_insert);
    return tc;
}

static TCase *
charz_tcase ()
{
    TCase *tc = tcase_create("char { 'z' }");
    tcase_add_checked_fixture(tc, charz_setup, teardown);
    tcase_add_test(tc, test_charz_init);
    tcase_add_test(tc, test_charz_append);
    tcase_add_test(tc, test_charz_prepend);
    tcase_add_test(tc, test_charz_insert);
    tcase_add_test(tc, test_charz_remove);
    return tc;
}

static TCase *
int246_tcase ()
{
    TCase *tc = tcase_create("int { 2, 4, 6 }");
    tcase_add_checked_fixture(tc, int246_setup, teardown);
    tcase_add_test(tc, test_int246_init);
    tcase_add_test(tc, test_int246_append);
    tcase_add_test(tc, test_int246_prepend);
    tcase_add_test(tc, test_int246_insert);
    tcase_add_test(tc, test_int246_remove);
    return tc;
}


Suite *
array_suite ()
{
    Suite *s = suite_create("array");
    suite_add_tcase(s, int_tcase());
    suite_add_tcase(s, charz_tcase());
    suite_add_tcase(s, int246_tcase());
    return s;
}
