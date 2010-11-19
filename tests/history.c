#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/history.h>

TCase *
history_tcase ()
{
    TCase *tc = tcase_create("history");
    return tc;
}

Suite *
history_suite ()
{
    Suite *s = suite_create("history");
    suite_add_tcase(s, history_tcase());
    return s;
}
