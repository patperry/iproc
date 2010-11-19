
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include <stdlib.h>
#include <check.h>

extern Suite *ieee754_suite ();
extern Suite *vector_suite ();
extern Suite *events_suite ();

int
main ()
{
    int nfail;
    SRunner *sr = srunner_create(NULL);
    srunner_add_suite(sr, ieee754_suite());
    srunner_add_suite(sr, vector_suite());
    srunner_add_suite(sr, events_suite());
    srunner_run_all(sr, CK_ENV);
    nfail = srunner_ntests_failed(sr);
    srunner_free(sr);
    return (nfail == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
}
