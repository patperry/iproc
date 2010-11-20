#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/array.h>


static void
setup ()
{
}

static void
teardown ()
{
}

START_TEST (test_init)
{
}
END_TEST

static TCase *
array_tcase ()
{
    TCase *tc = tcase_create("array");
    tcase_add_checked_fixture(tc, setup, teardown);
    tcase_add_test(tc, test_init);
    return tc;
}

Suite *
array_suite ()
{
    Suite *s = suite_create("array");
    suite_add_tcase(s, array_tcase());
    return s;
}
