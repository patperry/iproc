#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>
#include <setjmp.h>
#include <stdlib.h>

#include "../src/history.h"
#include "../src/design.h"
#include "../src/design2.h"

#include "coreutil.h"
#include "testutil.h"
#include "cmockery.h"
#include "ieee754.h"
#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"
#include "xalloc.h"

#include "fixtures/history.h"
#include "fixtures/design.h"
#include "fixtures/design2.h"


#include "fixtures/history.h"
#include "fixtures/design.h"


struct fixture {
	struct history_fixture h;
	struct design2_fixture d;
};

struct fixture ctx;

#define NRECV	ctx.h.nrecv
#define NSEND	ctx.h.nsend



struct test {
	struct history h;
	struct design2 d;
};

struct test test;

#define DESIGN2	&test.d
#define HISTORY	&test.h



static void fixture_setup_enron()
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	history_fixture_setup_enron(&ctx.h);
	design2_fixture_setup_enron(&ctx.d);
}

static void fixture_teardown()
{
	design2_fixture_teardown(&ctx.d);
	history_fixture_teardown(&ctx.h);
	print_message("\n\n");
}


static void setup()
{
	history_test_setup(&test.h, &ctx.h);
	design2_test_setup(&test.d, &test.h, &ctx.d);
}

static void teardown()
{
	design2_test_teardown(&test.d);
	history_test_teardown(&test.h);
}


static void test_isend1()
{
	double window = 1000;
	const struct var2 *v = design2_add_tvar(DESIGN2, "ISend", VAR2_ISEND, window);
	double t0 = 910930020;
	const double *x;
	const size_t *ind;
	size_t nz;

	assert_true(v);

	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz); /* -inf */
	assert_int_equal(nz, 0);

	history_advance(HISTORY, t0); /* t */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 0);

	history_advance(HISTORY, double_nextup(t0)); /* (t)+ */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 1);
	assert_int_equal(ind[0], 58);
	assert_real_identical(x[0], 1.0);

	design2_get_tvar_matrix(DESIGN2, 58, &x, &ind, &nz);
	assert_int_equal(nz, 0);

	history_advance(HISTORY, t0 + window); /* t + w */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 1);
	assert_int_equal(ind[0], 58);
	assert_real_identical(x[0], 1.0);

	design2_get_tvar_matrix(DESIGN2, 58, &x, &ind, &nz);
	assert_int_equal(nz, 0);

	history_advance(HISTORY, double_nextup(t0 + window)); /* (t + w)+ */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 1);
	assert_int_equal(ind[0], 58);
	assert_real_identical(x[0], 0.0);

	design2_get_tvar_matrix(DESIGN2, 58, &x, &ind, &nz);
	assert_int_equal(nz, 0);
}


static void test_isend2()
{
	double window = 2 * 7 * 24 * 60 * 60;
	const struct var2 *v = design2_add_tvar(DESIGN2, "ISend",
						VAR2_ISEND, window);
	size_t ito;
	size_t isend, jrecv;
	const struct message *msgs;
	size_t imsg, nmsg;

	history_get_messages(HISTORY, &msgs, &nmsg);


	double tlast[NSEND][NRECV];

	for (isend = 0; isend < NSEND; isend++) {
		for (jrecv = 0; jrecv < NRECV; jrecv++) {
			tlast[isend][jrecv] = -INFINITY;
		}
	}

	for (imsg = 0; imsg < nmsg; imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;
		history_advance(HISTORY, t);

		isend = msg->from;
		for (ito = 0; ito < msg->nto; ito++) {
			jrecv = msg->to[ito];
			const double *x = design2_tvar(DESIGN2, v, isend, jrecv);
			double tmsg = tlast[isend][jrecv];

			if (tmsg + window >= t) {
				assert_real_identical(x[0], 1.0);
			} else if (x) {
				assert_real_identical(x[0], 0.0);
			}
		}

		if (imsg + 1 < nmsg && msgs[imsg + 1].time != t) {
			while (msg->time == t) {
				for (ito = 0; ito < msg->nto; ito++) {
					tlast[msg->from][msg->to[ito]] = msg->time;
				}
				if (msg == msgs)
					break;
				msg--;
			}
		}
	}
}



static void test_irecv()
{
	double window = 5 * 7 * 24 * 60 * 60;
	const struct var2 *v = design2_add_tvar(DESIGN2, "IRecv",
						VAR2_IRECV, window);
	size_t ito;
	size_t isend, jrecv;
	const struct message *msgs;
	size_t imsg, nmsg;

	history_get_messages(HISTORY, &msgs, &nmsg);


	double tlast[NSEND][NRECV];

	for (isend = 0; isend < NSEND; isend++) {
		for (jrecv = 0; jrecv < NRECV; jrecv++) {
			tlast[isend][jrecv] = -INFINITY;
		}
	}

	for (imsg = 0; imsg < nmsg; imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;
		history_advance(HISTORY, t);

		jrecv = msg->from;
		for (ito = 0; ito < msg->nto; ito++) {
			isend = msg->to[ito];
			const double *x = design2_tvar(DESIGN2, v, isend, jrecv);
			double tmsg = tlast[jrecv][isend];

			if (tmsg + window >= t) {
				assert_real_identical(x[0], 1.0);
			} else if (x) {
				assert_real_identical(x[0], 0.0);
			}
		}

		if (imsg + 1 < nmsg && msgs[imsg + 1].time != t) {
			while (msg->time == t) {
				for (ito = 0; ito < msg->nto; ito++) {
					tlast[msg->from][msg->to[ito]] = msg->time;
				}
				if (msg == msgs)
					break;
				msg--;
			}
		}
	}
}


static void test_nsend()
{
	double intvls[] = { INFINITY };
	size_t nintvl = 1;
	const struct var2 *v = design2_add_tvar(DESIGN2, "NSend",
						VAR2_NSEND, intvls, nintvl);
	size_t ito;
	size_t isend, jrecv;
	const struct message *msgs;
	size_t imsg, nmsg;

	history_get_messages(HISTORY, &msgs, &nmsg);

	double xn[NSEND][NRECV];

	for (isend = 0; isend < NSEND; isend++) {
		for (jrecv = 0; jrecv < NRECV; jrecv++) {
			xn[isend][jrecv] = 0;
		}
	}

	for (imsg = 0; imsg < nmsg; imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;
		history_advance(HISTORY, t);

		isend = msg->from;
		for (ito = 0; ito < msg->nto; ito++) {
			jrecv = msg->to[ito];
			const double *x = design2_tvar(DESIGN2, v, isend, jrecv);
			double val = xn[isend][jrecv];

			if (x) {
				assert_real_identical(val, x[0]);
			} else {
				assert_real_identical(val, 0);
			}
		}

		if (imsg + 1 < nmsg && msgs[imsg + 1].time != t) {
			while (msg->time == t) {
				double wt = 1.0 / msg->nto;
				for (ito = 0; ito < msg->nto; ito++) {
					xn[msg->from][msg->to[ito]] += wt;
				}
				if (msg == msgs)
					break;
				msg--;
			}
		}
	}
}


static void test_nrecv()
{
	double intvls[] = { INFINITY };
	size_t nintvl = 1;
	const struct var2 *v = design2_add_tvar(DESIGN2, "NRecv",
						VAR2_NRECV, intvls, nintvl);
	size_t ito;
	size_t isend, jrecv;
	const struct message *msgs;
	size_t imsg, nmsg;

	history_get_messages(HISTORY, &msgs, &nmsg);

	double xn[NSEND][NRECV];

	for (isend = 0; isend < NSEND; isend++) {
		for (jrecv = 0; jrecv < NRECV; jrecv++) {
			xn[isend][jrecv] = 0;
		}
	}

	for (imsg = 0; imsg < nmsg; imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;
		history_advance(HISTORY, t);

		jrecv = msg->from;
		for (ito = 0; ito < msg->nto; ito++) {
			isend = msg->to[ito];
			const double *x = design2_tvar(DESIGN2, v, isend, jrecv);
			double val = xn[jrecv][isend];

			if (x) {
				assert_real_identical(val, x[0]);
			} else {
				assert_real_identical(val, 0);
			}
		}

		if (imsg + 1 < nmsg && msgs[imsg + 1].time != t) {
			while (msg->time == t) {
				double wt = 1.0 / msg->nto;
				for (ito = 0; ito < msg->nto; ito++) {
					xn[msg->from][msg->to[ito]] += wt;
				}
				if (msg == msgs)
					break;
				msg--;
			}
		}
	}
}


static void test_nsend2()
{
	double intvls1[] = { 4 * 60 * 60, 7 * 24 * 60 * 60, 4 * 7 * 24 * 60 * 60 };
	size_t nintvl1 = 3;
	double intvls2[] = { 2 * 60 * 60, 5 * 7 * 24 * 60 * 60 };
	size_t nintvl2 = 2;

	double *zero = xcalloc(nintvl1 * nintvl2, sizeof(double));
	const struct var2 *v1 = design2_add_tvar(DESIGN2, "NSend.1", VAR2_NSEND,
						 intvls1, nintvl1);
	const struct var2 *v2 = design2_add_tvar(DESIGN2, "NSend.2", VAR2_NSEND,
						 intvls2, nintvl2);
	const struct var2 *v = design2_add_tvar(DESIGN2, "NSend2", VAR2_NSEND2,
					       intvls1, nintvl1, intvls2, nintvl2);
	double *x = xcalloc(nintvl1 * nintvl2, sizeof(double));
	size_t i, j, h;

	assert_true(v);

	history_advance(HISTORY, 1003226742);

	size_t m = design2_count1(DESIGN2);
	size_t n = design2_count2(DESIGN2);
	assert(m == n);

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			memset(x, 0, nintvl1 * nintvl2 * sizeof(double));
			for (h = 0; h < n; h++) {
				const double *x1 = design2_tvar(DESIGN2, v1, i, h);
				const double *x2 = design2_tvar(DESIGN2, v2, h, j);
				x1 = design2_tvar(DESIGN2, v1, i, h); // second call to design2_tvar may reallocate

				if (x1 && x2) {
					blas_dger(nintvl2, nintvl1, 1.0, x2, 1, x1, 1, x, nintvl2);
				}
			}
			const double *xij = design2_tvar(DESIGN2, v, i, j);
			
			if (xij) {
				assert_vec_approx(xij, x, nintvl1 * nintvl2);
			} else {
				assert_vec_approx(zero, x, nintvl1 * nintvl2);
			}
		}
	}

	free(x);
	free(zero);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_isend1, setup, teardown),
		unit_test_setup_teardown(test_isend2, setup, teardown),
		unit_test_setup_teardown(test_irecv, setup, teardown),
		unit_test_setup_teardown(test_nsend, setup, teardown),
		unit_test_setup_teardown(test_nrecv, setup, teardown),
		unit_test_setup_teardown(test_nsend2, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown),

	};
	return run_tests(tests);
}


