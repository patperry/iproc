#include "port.h"
#include <math.h>
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "coreutil.h"
#include "cmockery.h"
#include "xalloc.h"

#include "enron.h"
#include "messages.h"
#include "design.h"
#include "var.h"
#include "frame.h"


static size_t nsend;
static size_t nrecv;
static size_t ntrait;
static double *traits;
static const char * const * trait_names;
static int has_loops;
static struct messages messages;
static struct design *recv_design;
static struct design2 *dyad_design;
static struct history *history;
static struct frame frame;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&nsend, &traits, &ntrait, &trait_names, ENRON_TERMS_MAX);
	nrecv = nsend;
	enron_messages_init(&messages, -1);
	has_loops = 0;	
}

static void enron_teardown_fixture()
{
	messages_deinit(&messages);
	free(traits);
	print_message("\n\n");
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
