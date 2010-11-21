#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/actors.h>

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
actors_tcase ()
{
    TCase *tc = tcase_create("actors");
    tcase_add_checked_fixture(tc, setup, teardown);
    tcase_add_test(tc, test_init);
    return tc;
}

Suite *
actors_suite ()
{
    Suite *s = suite_create("actors");
    suite_add_tcase(s, actors_tcase());
    return s;
}
