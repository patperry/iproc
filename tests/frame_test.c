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


static void dv_nsend_setup()
{
	double intvls[] = { INFINITY };
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 1);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, trait_names, traits, ntrait);

	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "NSend",VAR2_NSEND);

	history = frame_history(&frame);
}

static void dv_nsend_teardown()
{
	frame_deinit(&frame);
}

static void test_dv_nsend()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend, jrecv;
	const struct message *msg = NULL;
	struct messages_iter it;
	const struct var2 *v = design2_var(dyad_design, "NSend");
	
	double xnsend[nsend][nrecv];
	memset(xnsend, 0, nsend * nrecv * sizeof(double));

	isend = 0;

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		history_advance(history, t);

		isend = msg ? msg->from : 0;
		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			const double *dx = design2_tvar(dyad_design, v, isend, jrecv);
			
			if (dx) {
				assert(dx[0] == xnsend[isend][jrecv]);
				assert_true(dx[0] == xnsend[isend][jrecv]);
			} else {
				assert(0.0 == xnsend[isend][jrecv]);
				assert_true(0.0 == xnsend[isend][jrecv]);
			}
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			history_add(history, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				xnsend[msg->from][msg->to[ito]] += 1.0;
			}
		}
	}
}


static void dv_nrecv_setup()
{
	double intvls[] = { INFINITY };
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 1);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, trait_names, traits, ntrait);
	
	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "NRecv", VAR2_NRECV);

	history = frame_history(&frame);
}

static void dv_nrecv_teardown()
{
	frame_deinit(&frame);
}

static void test_dv_nrecv()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend;
	size_t jrecv, nrecv = frame_recv_count(&frame);
	const struct message *msg = NULL;
	struct messages_iter it;
	const struct var2 *v = design2_var(dyad_design, "NRecv");

	double xnsend[nsend][nrecv];
	memset(xnsend, 0, nsend * nrecv * sizeof(double));

	assert(nsend == nrecv);

	isend = 0;

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		history_advance(history, t);

		isend = msg ? msg->from : 0;

		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			const double *dx = design2_tvar(dyad_design, v, isend, jrecv);
			
			if (dx) {
				assert(dx[0] == xnsend[jrecv][isend]);
				assert_true(dx[0] == xnsend[jrecv][isend]);
			} else {
				assert_true(0.0 == xnsend[jrecv][isend]);
			}
		}

		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			history_add(history, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				xnsend[msg->from][msg->to[ito]] += 1.0;
			}
		}
	}
}

static void dv_irecv_setup()
{
	double intvls[4] = {
		112.50,  450.00, 1800.00, INFINITY
	};
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 4);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, trait_names, traits, ntrait);

	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "IRecv", VAR2_IRECV);

	history = frame_history(&frame);
}

static void dv_irecv_teardown()
{
	frame_deinit(&frame);
}


static void test_dv_irecv()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend;
	size_t jrecv, nrecv = frame_recv_count(&frame);
	const struct message *msg = NULL;
	struct messages_iter it;
	const struct var2 *v = design2_var(dyad_design, "IRecv");

	double tmsg;

	double tlast[nsend][nrecv];
	for (isend = 0; isend < nsend; isend++) {
		for (jrecv = 0; jrecv < nrecv; jrecv++) {
			tlast[isend][jrecv] = -INFINITY;
		}
	}

	assert(nsend == nrecv);


	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		history_advance(history, t);

		if (msg) {
			isend = msg->from;
			for (ito = 0; ito < msg->nto; ito++) {
				jrecv = msg->to[ito];
				const double *dx = design2_tvar(dyad_design, v, isend, jrecv);
				
				tmsg = tlast[jrecv][isend];
				
				if (isfinite(tmsg)) {
					assert(dx[0] == 1.0);
					assert_true(dx[0] == 1.0);
				} else {
					assert_true(!dx || dx[0] == 0.0);
				}
			}
		}

		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			history_add(history, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				tlast[msg->from][msg->to[ito]] = msg->time;
			}
		}
	}
}


static void dv_isend_setup()
{
	double intvls[4] = {
		112.50,  450.00, 1800.00, INFINITY
	};
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 4);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, trait_names, traits, ntrait);

	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "ISend", VAR2_ISEND);

	history = frame_history(&frame);
}


static void dv_isend_teardown()
{
	frame_deinit(&frame);
}


static void test_dv_isend()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend;
	size_t jrecv, nrecv = frame_recv_count(&frame);
	const struct message *msg = NULL;
	struct messages_iter it;
	const struct var2 *v = design2_var(dyad_design, "ISend");

	double tmsg;

	double tlast[nsend][nrecv];

	for (isend = 0; isend < nsend; isend++) {
		for (jrecv = 0; jrecv < nrecv; jrecv++) {
			tlast[isend][jrecv] = -INFINITY;
		}
	}

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		history_advance(history, t);

		if (msg) {
			isend = msg->from;
			for (ito = 0; ito < msg->nto; ito++) {
				jrecv = msg->to[ito];
				const double *dx = design2_tvar(dyad_design, v, isend, jrecv);
				
				tmsg = tlast[isend][jrecv];
				
				if (isfinite(tmsg)) {
					assert(dx[0] == 1.0);
					assert_true(dx[0] == 1.0);
				} else {
					assert_true(!dx || dx[0] == 0.0);
				}
			}
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			history_add(history, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				tlast[msg->from][msg->to[ito]] = msg->time;
			}
		}
	}
}




int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_dv_nsend, dv_nsend_setup, dv_nsend_teardown),
		unit_test_setup_teardown(test_dv_nrecv, dv_nrecv_setup, dv_nrecv_teardown),
		unit_test_setup_teardown(test_dv_irecv, dv_irecv_setup, dv_irecv_teardown),
		unit_test_setup_teardown(test_dv_isend, dv_isend_setup, dv_isend_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
