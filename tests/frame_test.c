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
static size_t ncohort;
static size_t ntrait;
static size_t *cohorts;
static struct dmatrix traits;
static const char * const * cohort_names;
static const char * const * trait_names;
static int has_effects;
static int has_loops;
static struct messages messages;
static struct design *design;
static size_t rv_irecv_index;
static size_t rv_isend_index;
static size_t rv_nrecv_index;
static size_t rv_nsend_index;
static struct frame frame;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&nsend,
			     &cohorts, &ncohort, &cohort_names,
			     &traits.data, &ntrait, &trait_names);
	traits.lda = MAX(1, nsend);
	nrecv = nsend;
	enron_messages_init(&messages, -1);
}

static void enron_teardown_fixture()
{
	free(cohorts);
	messages_deinit(&messages);
	free(traits.data);
	print_message("\n\n");
}

static void rv_nsend_setup()
{
	has_effects = 0;
	has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, NULL, 0);
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
	design_add_dvar(design, RECV_VAR_NSEND, NULL);
	rv_nsend_index = design_dvar_index(design, RECV_VAR_NSEND);
}

static void rv_nsend_teardown()
{
	frame_deinit(&frame);
}

static void test_rv_nsend()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend, jrecv;
	const struct message *msg = NULL;
	struct messages_iter it;
	double *x, *y;
	
	struct dmatrix xnsend = { xcalloc(nsend * nrecv, sizeof(double)), MAX(1, nsend) };
	
	x = xcalloc(design_dim(design), sizeof(double));
	x[rv_nsend_index] = 1.0;
	y = xmalloc(design_count(design) * sizeof(double));
	
	isend = 0;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);
		
		isend = msg ? msg->from : 0;
		frame_recv_mul(1.0, BLAS_NOTRANS, &frame, isend, x, 0.0, y);
		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			assert(y[jrecv] == MATRIX_ITEM(&xnsend, isend, jrecv));
			assert_true(y[jrecv] == MATRIX_ITEM(&xnsend, isend, jrecv));
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			
			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(&xnsend, msg->from, msg->to[ito]) += 1.0;
			}
		}
	}
	
	free(y);
	free(x);	
	free(xnsend.data);
}

static void rv_nrecv_setup()
{
	has_effects = 0;
	has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, NULL, 0);
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
	design_add_dvar(design, RECV_VAR_NRECV, NULL);
	rv_nrecv_index = design_dvar_index(design, RECV_VAR_NRECV);
}

static void rv_nrecv_teardown()
{
	frame_deinit(&frame);
}

static void test_rv_nrecv()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend;
	size_t jrecv, nrecv = design_count(design);
	const struct message *msg = NULL;
	struct messages_iter it;

	double *x, *y;
	
	struct dmatrix xnrecv = { xcalloc(nsend * design_count(design), sizeof(double)),
				  MAX(1, nsend) };

	x = xcalloc(design_dim(design), sizeof(double));
	x[rv_nrecv_index] = 1.0;
	y = xmalloc(design_count(design) * sizeof(double));
	
	isend = 0;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);
		
		isend = msg ? msg->from : 0;
		frame_recv_mul(1.0, BLAS_NOTRANS, &frame, isend, x, 0.0, y);
		for (jrecv = 0; jrecv < nrecv; jrecv += 5) {
			assert(y[jrecv] == MATRIX_ITEM(&xnrecv, jrecv, isend));
			assert_true(y[jrecv] == MATRIX_ITEM(&xnrecv, jrecv, isend));
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			
			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(&xnrecv, msg->from, msg->to[ito]) += 1.0;
			}
		}
	}
	
	free(y);
	free(x);	
	free(xnrecv.data);
}

static void rv_irecv_setup()
{
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	has_effects = 0;
	has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
	design_add_dvar(design, RECV_VAR_IRECV, NULL);
	rv_irecv_index = design_dvar_index(design, RECV_VAR_IRECV);
}

static void rv_irecv_teardown()
{
	frame_deinit(&frame);
}


static void test_rv_irecv()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend;
	size_t jrecv, j, nrecv = design_count(design);
	const struct message *msg = NULL;
	struct messages_iter it;

	double x = 1.0;
	size_t k = 0;
	struct vpattern pat_k;
	pat_k.indx = &k;
	pat_k.nz = 1;
	double *y;
	double tmsg;
	
	struct dmatrix tlast = { xmalloc(nsend * design_count(design) * sizeof(double)),
		MAX(1, nsend) };
	for (j = 0; j < nsend * design_count(design); j++) {
		tlast.data[j] = -INFINITY;
	}
	
	y = xmalloc(design_count(design) * sizeof(double));
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);
		
		isend = msg ? msg->from : 0;
		jrecv = msg ? msg->to[0] : 0;
		
		k = 0;
		frame_recv_dmuls(1.0, BLAS_NOTRANS, &frame, isend, &x, &pat_k, 0.0, y);

		for (j = 0; j < 5; j++) {
			size_t ix = (jrecv + j) % nrecv;
			tmsg = MATRIX_ITEM(&tlast, ix, isend);
			if (isfinite(tmsg)) {
				assert(y[ix] == 1.0);
				assert_true(y[ix] == 1.0);
			} else {
				assert_true(y[ix] == 0.0);
			}
		}
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			
			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(&tlast, msg->from, msg->to[ito]) = msg->time;
			}
		}
	}

	free(y);	
	free(tlast.data);
}


static void rv_isend_setup()
{
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	has_effects = 0;
	has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
	design_add_dvar(design, RECV_VAR_ISEND, NULL);
	rv_isend_index = design_dvar_index(design, RECV_VAR_ISEND);
}

static void rv_isend_teardown()
{
	frame_deinit(&frame);
}


static void test_rv_isend()
{
	double t;
	size_t itie, ntie, ito;
	size_t isend;
	size_t jrecv, j, nrecv = design_count(design);
	const struct message *msg = NULL;
	struct messages_iter it;

	double x = 1.0;
	size_t k = 0;
	struct vpattern pat_k;
	pat_k.indx = &k;
	pat_k.nz = 1;
	double *y;
	double tmsg;

	struct dmatrix tlast = { xmalloc(nsend * design_count(design) * sizeof(double)), MAX(1, nsend) };
	for (j = 0; j < nsend * design_count(design); j++) {
		tlast.data[j] = -INFINITY;
	}

	y = xmalloc(design_count(design) * sizeof(double));

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		frame_advance(&frame, t);

		isend = msg ? msg->from : 0;
		jrecv = msg ? msg->to[0] : 0;
		k = 0;
		frame_recv_dmuls(1.0, BLAS_NOTRANS, &frame, isend, &x, &pat_k, 0.0, y);

		for (j = 0; j < 5; j++) {
			size_t ix = (jrecv + j) % nrecv;
			tmsg = MATRIX_ITEM(&tlast, isend, ix);
			if (isfinite(tmsg)) {
				assert(y[ix] == 1.0);
				assert_true(y[ix] == 1.0);
			} else {
				assert_true(y[ix] == 0.0);
			}
		}
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);

			for (ito = 0; ito < msg->nto; ito++) {
				MATRIX_ITEM(&tlast, msg->from, msg->to[ito]) = msg->time;
			}
		}
	}

	free(y);
	free(tlast.data);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_rv_nsend, rv_nsend_setup, rv_nsend_teardown),
		unit_test_setup_teardown(test_rv_nrecv, rv_nrecv_setup, rv_nrecv_teardown),
		unit_test_setup_teardown(test_rv_irecv, rv_irecv_setup, rv_irecv_teardown),
		unit_test_setup_teardown(test_rv_isend, rv_isend_setup, rv_isend_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
