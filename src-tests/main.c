
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdlib.h>
#include <check.h>

extern Suite *actors_suite ();
extern Suite *array_suite ();
extern Suite *ieee754_suite ();
extern Suite *vector_suite ();
extern Suite *svector_suite ();
extern Suite *events_suite ();
extern Suite *history_suite ();

int
main ()
{
    int nfail;
    SRunner *sr = srunner_create(NULL);
    srunner_add_suite(sr, array_suite());
    srunner_add_suite(sr, ieee754_suite());
    srunner_add_suite(sr, vector_suite());
    srunner_add_suite(sr, svector_suite());
    srunner_add_suite(sr, events_suite());
    srunner_add_suite(sr, history_suite());
    srunner_add_suite(sr, actors_suite());
    srunner_run_all(sr, CK_ENV);
    nfail = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (nfail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}