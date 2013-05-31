#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"

#include "ieee754.h"
#include "enron/messages.h"
#include "testutil.h"
#include "history.h"

static size_t NSEND;
static size_t NRECV;
static double *TIME;
static size_t *FROM;
static size_t **TO;
static size_t *NTO;
static intptr_t *ATTR;
static size_t NMSG;

static struct history HISTORY;


static void enron_setup_fixture()
{
	int ok = enron_messages_init(&NSEND, &NRECV, &TIME, &FROM, &TO, &NTO,
				     &ATTR, &NMSG, 0);
	assert(ok);


}

static void enron_teardown_fixture()
{
	size_t i;

	for (i = 0; i < NMSG; i++) {
		free(TO[i]);
	}

	free(ATTR);
	free(NTO);
	free(TO);
	free(FROM);
	free(TIME);
}

static void enron_setup()
{
	size_t i, n = NMSG;
	struct history *h = &HISTORY;

	history_init(h, NSEND, NRECV);
	for (i = 0; i < n; i++) {
		history_advance(h, TIME[i]);
		history_add(h, FROM[i], TO[i], NTO[i], ATTR[i]);
	}
	history_reset(h);
}

static void enron_teardown()
{
	history_deinit(&HISTORY);
}

static void test_sizes()
{
	assert_int_equal(history_nsend(&HISTORY), NSEND);
	assert_int_equal(history_nrecv(&HISTORY), NRECV);
}

static void test_messages()
{
	struct message *msg;
	size_t i;

	history_advance(&HISTORY, INFINITY);

	assert_int_equal(history_count(&HISTORY), NMSG);
	for (i = 0; i < NMSG; i++) {
		msg = history_item(&HISTORY, i);
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

	history_advance(&HISTORY, t);

	assert_real_identical(history_time(&HISTORY), t);
	assert_int_equal(history_count(&HISTORY), 0);

	for (i = 0; i < NSEND; i++) {
		struct history_actor *ha = history_send(&HISTORY, i);
		history_actor_get_msgs(ha, &ind, &len);
		assert_int_equal(len, 0);

		history_actor_get_alters(ha, &wts, &ind, &len);
		assert_int_equal(len, 0);
	}

	for (j = 0; j < NRECV; j++) {
		struct history_actor *ha = history_recv(&HISTORY, j);
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

	history_advance(&HISTORY, t);

	assert_real_identical(history_time(&HISTORY), t);
	assert_int_equal(history_count(&HISTORY), 1);

	for (i = 0; i < NSEND; i++) {
		struct history_actor *ha = history_send(&HISTORY, i);
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
		struct history_actor *ha = history_recv(&HISTORY, j);
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
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_sizes, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_messages, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_advance_tfirst, enron_setup, enron_teardown),
		unit_test_setup_teardown(test_advance_after_tfirst, enron_setup, enron_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
