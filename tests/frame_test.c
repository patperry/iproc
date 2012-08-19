#include "port.h"
#include <math.h>
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "coreutil.h"
#include "cmockery.h"
#include "xalloc.h"

#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"


static size_t nsend;
static size_t nrecv;
static size_t ntrait;
static double *traits;
static size_t ldtraits;
static const char * const * trait_names;
static int has_loops;
static struct messages messages;
static struct design *recv_design;
static struct design2 *dyad_design;
static struct frame frame;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&nsend, &traits, &ntrait, &trait_names);
	ldtraits = MAX(1, nsend);
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
	frame_init(&frame, nsend, nrecv, has_loops, NULL, 0);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, ntrait, trait_names, traits);

	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "NSend", DYAD_VAR_NSEND);
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
	
	double *xnsend = xcalloc(nsend * nrecv, sizeof(double));
	size_t ldxnsend = MAX(1, nsend);

	isend = 0;

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);

		isend = msg ? msg->from : 0;
		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			const double *dx = design2_tvars(dyad_design, isend, jrecv);
			
			if (dx) {
				assert(dx[0] == MATRIX_ITEM(xnsend, ldxnsend, isend, jrecv));
				assert_true(dx[0] == MATRIX_ITEM(xnsend, ldxnsend, isend, jrecv));
			} else {
				assert_true(0.0 == MATRIX_ITEM(xnsend, ldxnsend, isend, jrecv));
			}
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(xnsend, ldxnsend, msg->from, msg->to[ito]) += 1.0;
			}
		}
	}

	free(xnsend);
}


static void dv_nrecv_setup()
{
	frame_init(&frame, nsend, nrecv, has_loops, NULL, 0);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, ntrait, trait_names, traits);
	
	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "NRecv", DYAD_VAR_NRECV);
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

	double *xnsend = xcalloc(nsend * nrecv, sizeof(double));
	size_t ldxnsend = MAX(1, nsend);

	isend = 0;

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);

		isend = msg ? msg->from : 0;

		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			const double *dx = design2_tvars(dyad_design, isend, jrecv);
			
			if (dx) {
				assert(dx[0] == MATRIX_ITEM(xnsend, ldxnsend, jrecv, isend));
				assert_true(dx[0] == MATRIX_ITEM(xnsend, ldxnsend, jrecv, isend));
			} else {
				assert_true(0.0 == MATRIX_ITEM(xnsend, ldxnsend, jrecv, isend));
			}
		}

		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(xnsend, ldxnsend, msg->from, msg->to[ito]) += 1.0;
			}
		}
	}

	free(xnsend);
}

static void dv_irecv_setup()
{
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, ntrait, trait_names, traits);

	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "IRecv", DYAD_VAR_IRECV);
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
	size_t jrecv, j, nrecv = frame_recv_count(&frame);
	const struct message *msg = NULL;
	struct messages_iter it;

	double tmsg;

	double *tlast = xmalloc(nsend * nrecv * sizeof(double));
	size_t ldtlast = MAX(1, nsend);
	for (j = 0; j < nsend * nrecv; j++) {
		tlast[j] = -INFINITY;
	}

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);

		if (msg) {
			isend = msg->from;
			for (ito = 0; ito < msg->nto; ito++) {
				jrecv = msg->to[ito];
				const double *dx = design2_tvars(dyad_design, isend, jrecv);
				
				tmsg = MATRIX_ITEM(tlast, ldtlast, jrecv % nrecv, isend);
				
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
			frame_add(&frame, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(tlast, ldtlast, msg->from, msg->to[ito]) = msg->time;
			}
		}
	}

	free(tlast);
}


static void dv_isend_setup()
{
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);

	recv_design = frame_recv_design(&frame);
	design_add_traits(recv_design, ntrait, trait_names, traits);

	dyad_design = frame_dyad_design(&frame);
	design2_add_tvar(dyad_design, "ISend", DYAD_VAR_ISEND);
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
	size_t jrecv, j, nrecv = frame_recv_count(&frame);
	const struct message *msg = NULL;
	struct messages_iter it;

	double tmsg;

	double *tlast = xmalloc(nsend * nrecv * sizeof(double));
	size_t ldtlast = MAX(1, nsend);
	for (j = 0; j < nsend * nrecv; j++) {
		tlast[j] = -INFINITY;
	}

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);

		if (msg) {
			isend = msg->from;
			for (ito = 0; ito < msg->nto; ito++) {
				jrecv = msg->to[ito];
				const double *dx = design2_tvars(dyad_design, isend, jrecv);
				
				tmsg = MATRIX_ITEM(tlast, ldtlast, isend, jrecv % nrecv);
				
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
			frame_add(&frame, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(tlast, ldtlast, msg->from, msg->to[ito]) = msg->time;
			}
		}
	}

	free(tlast);
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
