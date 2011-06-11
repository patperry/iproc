#include "port.h"
#include <math.h>
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
#include "vrecv.h"
#include "frame.h"


static struct actors senders;
static struct actors receivers;
static bool has_reffects;
static bool has_loops;
static struct vector intervals;
static struct messages messages;
static struct design design;
static ssize_t vnrecv_index;
static ssize_t vrecv_index;
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

static void vnrecv_setup(void **state)
{
	has_reffects = false;
	has_loops = false;
	vector_init(&intervals, 0);
	design_init(&design, &senders, &receivers, &intervals);
	design_set_loops(&design, has_loops);
	design_set_reffects(&design, has_reffects);
	design_add_var(&design, VAR_TYPE_NRECV);
	vnrecv_index = design_var_index(&design, VAR_TYPE_NRECV);
	frame_init(&frame, &design);
}

static void vnrecv_teardown(void **state)
{
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}



static void test_vnrecv(void **state)
{
	double t;
	ssize_t itie, ntie, ito;
	ssize_t isend;
	ssize_t jrecv, nrecv = design_nreceiver(&design);
	const struct message *msg = NULL;
	struct messages_iter it;
	struct matrix xnrecv;
	struct vector x, y;
	
	
	
	matrix_init(&xnrecv, design_nsender(&design), design_nreceiver(&design));
	matrix_fill(&xnrecv, 0.0);

	vector_init(&x, design_dim(&design));
	vector_set_basis(&x, vnrecv_index);
	vector_init(&y, design_nreceiver(&design));
	
	isend = 0;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance_to(&frame, t);
		
		isend = msg ? msg->from : 0;
		frame_mul(1.0, TRANS_NOTRANS, &frame, isend, &x, 0.0, &y);
		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			assert_true(vector_item(&y, jrecv) == matrix_item(&xnrecv, jrecv, isend));
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_insert(&frame, msg);
			
			for (ito = 0; ito < msg->nto; ito++) {
				*matrix_item_ptr(&xnrecv, msg->from, msg->to[ito]) += 1.0;
			}
		}
	}
	
	matrix_deinit(&xnrecv);
}


static void vrecv_setup(void **state)
{
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	has_reffects = false;
	has_loops = false;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	design_init(&design, &senders, &receivers, &intervals);
	design_set_loops(&design, has_loops);
	design_set_reffects(&design, has_reffects);
	design_add_var(&design, VAR_TYPE_RECV);
	vrecv_index = design_var_index(&design, VAR_TYPE_RECV);
	frame_init(&frame, &design);
}

static void vrecv_teardown(void **state)
{
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}



static void test_vrecv(void **state)
{
	double t;
	ssize_t itie, ntie, ito;
	ssize_t isend;
	ssize_t jrecv, j, nrecv = design_nreceiver(&design);
	const struct message *msg = NULL;
	struct messages_iter it;
	struct matrix tlast;
	struct svector x, y;
	double delta, tmsg, tlo, thi;
	ssize_t i, n = vector_dim(&intervals);
	
	matrix_init(&tlast, design_nsender(&design), design_nreceiver(&design));
	matrix_fill(&tlast, -INFINITY);
	
	svector_init(&x, design_dim(&design));
	svector_init(&y, design_nreceiver(&design));
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance_to(&frame, t);
		
		isend = msg ? msg->from : 0;
		jrecv = msg ? msg->to[0] : 0;
		
		for (i = 0; i <= n; i++) {
			tlo = i == 0 ? 0 : vector_item(&intervals, i - 1);
			thi = i == n ? INFINITY : vector_item(&intervals, i);
			
			svector_set_basis(&x, vrecv_index + i);
			frame_dmuls(1.0, TRANS_NOTRANS, &frame, isend, &x, 0.0, &y);

			for (j = 0; j < 1; j++) {
				tmsg = matrix_item(&tlast, (jrecv + j) % nrecv, isend);
				delta = t - tmsg;
				if (isfinite(tmsg) && tlo < delta && delta <= thi) {
					assert_true(svector_item(&y, jrecv) == 1.0);
				} else {
					assert_true(svector_item(&y, (jrecv + j) % nrecv) == 0.0);
				}
			}
		}
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_insert(&frame, msg);
			
			for (ito = 0; ito < msg->nto; ito++) {
				*matrix_item_ptr(&tlast, msg->from, msg->to[ito]) = msg->time;
			}
		}
	}
	
	matrix_deinit(&tlast);
}


int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_vnrecv, vnrecv_setup, vnrecv_teardown),
		unit_test_setup_teardown(test_vrecv, vrecv_setup, vrecv_teardown),		
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
