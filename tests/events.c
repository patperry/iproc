#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/events.h>


Suite *
events_suite ()
{
    Suite *s = suite_create("events");
    return s;
}
