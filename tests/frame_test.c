#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include "cmockery.h"

#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vnrecv.h"
#include "frame.h"


static struct actors senders;
static struct actors receivers;
static bool has_reffects;
static struct messages messages;
static struct design design;
static struct vnrecv vnrecv;
static struct frame frame;


static void enron_setup_fixture(void **state)
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&senders);
	enron_employees_init(&receivers);
	enron_messages_init(&messages);
}

static void enron_teardown_fixture(void **state)
{
	messages_deinit(&messages);
	actors_deinit(&receivers);
	actors_deinit(&senders);
	print_message("\n\n");
}

static void enron_setup(void **state)
{
	has_reffects = false;
	design_init(&design, &senders, &receivers, has_reffects);
	vnrecv_init(&vnrecv, NULL, 0);
	frame_init(&frame, &design);
}

static void enron_teardown(void **state)
{
	frame_deinit(&frame);
	vnrecv_deinit(&vnrecv);
	design_deinit(&design);
}

static void test_basic(void **state)
{
	double t;
	ssize_t itie, ntie;
	const struct message *msg;
	struct messages_iter it = messages_iter(&messages);	
	
	while (messages_iter_advance(&it)) {
		t = messages_iter_current_time(&it);
		frame_advance_to(&frame, t, NULL);
		
		ntie = messages_iter_ntie(&it);
		for (itie = 0; itie < ntie; itie++) {
			msg = messages_iter_current(&it, itie);
			frame_insert(&frame, msg);
		}
	}
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_basic, enron_setup, enron_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
