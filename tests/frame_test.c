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
static bool has_loops;
static struct messages messages;
static struct design design;
static struct vnrecv vnrecv;
static ssize_t vnrecv_index;
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
	has_loops = false;
	design_init(&design, &senders, &receivers, has_reffects, has_loops);
	vnrecv_init(&vnrecv, NULL, 0);
	vnrecv_index = design_add_dyad_var(&design, &vnrecv.dyad_var);	
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
	ssize_t itie, ntie, ito;
	ssize_t isend;
	ssize_t jrecv, nrecv = design_nreceiver(&design);
	const struct message *msg = NULL;
	struct messages_iter it = messages_iter(&messages);	
	struct matrix xnrecv;
	struct vector x, y;
	
	
	
	matrix_init(&xnrecv, design_nsender(&design), design_nreceiver(&design));
	matrix_fill(&xnrecv, 0.0);

	vector_init(&x, design_dim(&design));
	vector_set_basis(&x, vnrecv_index);
	vector_init(&y, design_nreceiver(&design));
	
	isend = 0;
	
	while (messages_iter_advance(&it)) {
			//print_message(".");
			//fflush(stdout);
		t = messages_iter_current_time(&it);
		frame_advance_to(&frame, t, NULL);
		
		isend = msg ? msg->from : 0;
		frame_mul(1.0, TRANS_NOTRANS, &frame, isend, &x, 0.0, &y);
		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			assert_true(vector_get(&y, jrecv) == matrix_get(&xnrecv, jrecv, isend));
		}
		
		ntie = messages_iter_ntie(&it);
		for (itie = 0; itie < ntie; itie++) {
			msg = messages_iter_current(&it, itie);
			frame_insert(&frame, msg);
			
			for (ito = 0; ito < msg->nto; ito++) {
				*matrix_at(&xnrecv, msg->from, msg->to[ito]) += 1.0;
			}
		}
	}
	
	matrix_deinit(&xnrecv);
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
