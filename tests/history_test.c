#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>

#include <history.h>

#include "cmockery.h"
#include "ieee754.h"
#include "testutil.h"

#include "fixtures/history.h"


struct fixture {
	struct history_fixture h;
};

struct fixture ctx;

#define NSEND	ctx.h.nsend
#define NRECV	ctx.h.nrecv
#define TIME	ctx.h.time
#define FROM	ctx.h.from
#define TO	ctx.h.to
#define NTO	ctx.h.nto
#define ATTR	ctx.h.attr
#define NMSG	ctx.h.count


struct test {
	struct history h;
};

struct test test;

#define HISTORY	&test.h


static void fixture_setup_enron()
{
	history_fixture_setup_enron(&ctx.h);
}

static void fixture_teardown()
{
	history_fixture_teardown(&ctx.h);
}

static void setup()
{
	history_test_setup(&test.h, &ctx.h);
}

static void teardown()
{
	history_test_teardown(&test.h);
}

static void test_sizes()
{
	assert_int_equal(history_nsend(HISTORY), NSEND);
	assert_int_equal(history_nrecv(HISTORY), NRECV);
}

static void test_messages()
{
	struct message *msg;
	size_t i;

	history_advance(HISTORY, INFINITY);

	assert_int_equal(history_count(HISTORY), NMSG);
	for (i = 0; i < NMSG; i++) {
		msg = history_item(HISTORY, i);
		assert_real_identical(msg->time, TIME[i]);
		assert_int_equal(msg->from, FROM[i]);
		assert_int_equal(msg->nto, NTO[i]);
		assert_memory_equal(msg->to, TO[i], NTO[i] * sizeof(size_t));
		assert_int_equal(msg->attr, ATTR[i]);
	}

}

static void test_advance_tfirst()
{
	double t = TIME[0];
	const size_t *ind;
	const double *wts;
	size_t i, j, len;

	history_advance(HISTORY, t);

	assert_real_identical(history_time(HISTORY), t);
	assert_int_equal(history_count(HISTORY), 0);

	for (i = 0; i < NSEND; i++) {
		struct history_actor *ha = history_send(HISTORY, i);
		history_actor_get_msgs(ha, &ind, &len);
		assert_int_equal(len, 0);

		history_actor_get_alters(ha, &wts, &ind, &len);
		assert_int_equal(len, 0);
	}

	for (j = 0; j < NRECV; j++) {
		struct history_actor *ha = history_recv(HISTORY, j);
		history_actor_get_msgs(ha, &ind, &len);
		assert_int_equal(len, 0);

		history_actor_get_alters(ha, &wts, &ind, &len);
		assert_int_equal(len, 0);
	}
}

static void test_advance_after_tfirst()
{
	double t = double_nextup(TIME[0]);
	const size_t *ind, *mind;
	const double *wts;
	size_t i, j, len, mlen;

	assert(TIME[0] != TIME[1]);
	assert(NTO[0] == 1);

	history_advance(HISTORY, t);

	assert_real_identical(history_time(HISTORY), t);
	assert_int_equal(history_count(HISTORY), 1);

	for (i = 0; i < NSEND; i++) {
		struct history_actor *ha = history_send(HISTORY, i);
		history_actor_get_msgs(ha, &mind, &mlen);
		history_actor_get_alters(ha, &wts, &ind, &len);

		if (i != FROM[0]) {
			assert_int_equal(mlen, 0);
			assert_int_equal(mlen, 0);
		} else {
			assert_int_equal(mlen, 1);
			assert_int_equal(mind[0], 0);

			assert_int_equal(len, NTO[0]);
			assert_real_identical(wts[0], 1.0);
			assert_int_equal(ind[0], TO[0][0]);
		}
	}

	for (j = 0; j < NRECV; j++) {
		struct history_actor *ha = history_recv(HISTORY, j);
		history_actor_get_msgs(ha, &mind, &mlen);
		history_actor_get_alters(ha, &wts, &ind, &len);

		if (j != TO[0][0]) {
			assert_int_equal(mlen, 0);
			assert_int_equal(mlen, 0);
		} else {
			assert_int_equal(mlen, 1);
			assert_int_equal(mind[0], 0);

			assert_int_equal(len, 1);
			assert_real_identical(wts[0], 1.0 / NTO[0]);
			assert_int_equal(ind[0], FROM[0]);
		}

	}
}


int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_sizes, setup, teardown),
		unit_test_setup_teardown(test_messages, setup, teardown),
		unit_test_setup_teardown(test_advance_tfirst, setup, teardown),
		unit_test_setup_teardown(test_advance_after_tfirst, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown),
	};
	return run_tests(tests);
}
