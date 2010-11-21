#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/history.h>

iproc_history *history;

static void
setup ()
{
    history = iproc_history_new();
}

static void
teardown ()
{
    iproc_history_free(history);
}

START_TEST (test_init)
{
    ck_assert(history != NULL);
    ck_assert_int_eq(iproc_history_nsend(history), 0);
    ck_assert_int_eq(iproc_history_nrecv(history), 0);
}
END_TEST

START_TEST (test_insert_empty)
{
    iproc_history_insert(history, 7, 3);
    ck_assert_int_eq(iproc_events_ncur(iproc_history_send(history, 0)), 0);
    ck_assert_int_eq(iproc_events_ncur(iproc_history_send(history, 6)), 0);
    ck_assert_int_eq(iproc_events_ncur(iproc_history_recv(history, 0)), 0);
    ck_assert_int_eq(iproc_events_ncur(iproc_history_recv(history, 2)), 0);
}
END_TEST

START_TEST (test_insert1)
{
    iproc_history_insert(history, 2, 10);
    ck_assert_int_eq(iproc_history_nsend(history), 3);
    ck_assert_int_eq(iproc_history_nrecv(history), 11);
    ck_assert(iproc_events_cur(iproc_history_send(history, 2), 0) == 10);
    ck_assert(iproc_events_cur(iproc_history_recv(history, 10), 0) == 2);
}
END_TEST

START_TEST (test_insert2)
{
    iproc_history_insert(history, 3, 1);
    iproc_history_insert(history, 2, 4);
    ck_assert_int_eq(iproc_history_nsend(history), 4);
    ck_assert_int_eq(iproc_history_nrecv(history), 5);

    iproc_history_insert(history, 3, 2);
    ck_assert_int_eq(iproc_history_nsend(history), 4);
    ck_assert_int_eq(iproc_history_nrecv(history), 5);

    ck_assert(iproc_events_cur(iproc_history_send(history, 2), 0) == 4);
    ck_assert(iproc_events_cur(iproc_history_send(history, 3), 0) == 1);
    ck_assert(iproc_events_cur(iproc_history_send(history, 3), 1) == 2);

    ck_assert(iproc_events_cur(iproc_history_recv(history, 1), 0) == 3);
    ck_assert(iproc_events_cur(iproc_history_recv(history, 2), 0) == 3);
    ck_assert(iproc_events_cur(iproc_history_recv(history, 4), 0) == 2);
}
END_TEST

START_TEST (test_insert_advance)
{
    iproc_history_insert(history, 8, 0);
    iproc_history_advance(history, 123);

    ck_assert_int_eq(iproc_events_past(iproc_history_send(history, 8), 0), 0);
    ck_assert_int_eq(iproc_events_past_dt(iproc_history_send(history, 8), 0), 123);
    ck_assert_int_eq(iproc_events_past(iproc_history_recv(history, 0), 0), 8);
    ck_assert_int_eq(iproc_events_past_dt(iproc_history_recv(history, 0), 0), 123);
}
END_TEST

START_TEST (test_clear)
{
    iproc_history_insert(history, 0, 5);
    iproc_history_advance(history, 99);
    iproc_history_send(history, 0);
    iproc_history_clear(history);
    ck_assert_int_eq(iproc_history_nsend(history), 1);
    ck_assert_int_eq(iproc_history_nrecv(history), 6);

    iproc_history_insert(history, 0, 5);
    ck_assert_int_eq(iproc_events_cur(iproc_history_send(history, 0), 0), 5);
    ck_assert_int_eq(iproc_events_npast(iproc_history_send(history, 0)), 0);

    iproc_history_advance(history, 1);
    ck_assert_int_eq(iproc_events_past_dt(iproc_history_send(history, 0), 0), 1);
}
END_TEST

START_TEST (test_insertm)
{
    int64_t from[2] = { 3, 0 };
    int64_t to[3]   = { 2, 0, 3 };

    iproc_history_insertm(history, 2, from, 3, to);
    iproc_history_advance(history, 99);

    iproc_events *send0 = iproc_history_send(history, 0);
    iproc_events *send3 = iproc_history_send(history, 3);
    iproc_events *recv0 = iproc_history_recv(history, 0);
    iproc_events *recv2 = iproc_history_recv(history, 2);
    iproc_events *recv3 = iproc_history_recv(history, 3);

    ck_assert_int_eq(iproc_events_past(send0, 0), 0);
    ck_assert_int_eq(iproc_events_past(send0, 1), 2);
    ck_assert_int_eq(iproc_events_past(send0, 2), 3);

    ck_assert_int_eq(iproc_events_past(send3, 0), 0);
    ck_assert_int_eq(iproc_events_past(send3, 1), 2);
    ck_assert_int_eq(iproc_events_past(send3, 2), 3);

    ck_assert_int_eq(iproc_events_past(recv0, 0), 0);
    ck_assert_int_eq(iproc_events_past(recv0, 1), 3);

    ck_assert_int_eq(iproc_events_past(recv2, 0), 0);
    ck_assert_int_eq(iproc_events_past(recv2, 1), 3);

    ck_assert_int_eq(iproc_events_past(recv3, 0), 0);
    ck_assert_int_eq(iproc_events_past(recv3, 1), 3);
}
END_TEST


TCase *
history_tcase ()
{
    TCase *tc = tcase_create("history");
    tcase_add_checked_fixture(tc, setup, teardown);
    tcase_add_test(tc, test_init);
    tcase_add_test(tc, test_insert_empty);
    tcase_add_test(tc, test_insert1);
    tcase_add_test(tc, test_insert2);
    tcase_add_test(tc, test_insert_advance);
    tcase_add_test(tc, test_insertm);
    tcase_add_test(tc, test_clear);
    return tc;
}

Suite *
history_suite ()
{
    Suite *s = suite_create("history");
    suite_add_tcase(s, history_tcase());
    return s;
}
