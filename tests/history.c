#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/history.h>

iproc_history *history;

static void
setup ()
{
    history = iproc_history_new();
}

static void
teardown ()
{
    iproc_history_free(history);
}

TCase *
history_tcase ()
{
    TCase *tc = tcase_create("history");
    tcase_add_checked_fixture(tc, setup, teardown);
    return tc;
}

Suite *
history_suite ()
{
    Suite *s = suite_create("history");
    suite_add_tcase(s, history_tcase());
    return s;
}
