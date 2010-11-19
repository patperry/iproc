#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <check.h>
#include <iproc/events.h>


iproc_events *events;

void
setup ()
{
    events = iproc_events_new();
}

void
teardown ()
{
    iproc_events_free(events);
}

START_TEST (test_init)
{
    ck_assert_int_eq(iproc_events_ncur(events), 0);
    ck_assert_int_eq(iproc_events_npast(events), 0);
    ck_assert(iproc_events_find_cur(events, 0) < 0);
    ck_assert(iproc_events_find_past(events, 0) < 0);
}
END_TEST

START_TEST (test_insert1)
{
    iproc_events_insert(events, 12);
    int64_t i = iproc_events_find_cur(events, 12);

    ck_assert_int_eq(iproc_events_ncur(events), 1);
    ck_assert(i >= 0);
    ck_assert_int_eq(iproc_events_cur(events, i), 12);
    ck_assert_int_eq(iproc_events_npast(events), 0);
}
END_TEST

START_TEST (test_insert2)
{
    iproc_events_insert(events, 4);
    iproc_events_insert(events, 19);
    int64_t i1 = iproc_events_find_cur(events, 4);
    int64_t i2 = iproc_events_find_cur(events, 19);

    ck_assert_int_eq(iproc_events_ncur(events), 2);
    ck_assert_int_eq(iproc_events_cur(events, i1), 4);
    ck_assert_int_eq(iproc_events_cur(events, i2), 19);
    ck_assert_int_eq(iproc_events_npast(events), 0);

}
END_TEST

START_TEST (test_insert2_dup)
{
    iproc_events_insert(events, 4);
    iproc_events_insert(events, 19);
    iproc_events_insert(events, 4);
    iproc_events_insert(events, 19);
    int64_t i1 = iproc_events_find_cur(events, 4);
    int64_t i2 = iproc_events_find_cur(events, 19);

    ck_assert_int_eq(iproc_events_ncur(events), 2);
    ck_assert_int_eq(iproc_events_cur(events, i1), 4);
    ck_assert_int_eq(iproc_events_cur(events, i2), 19);
    ck_assert_int_eq(iproc_events_npast(events), 0);
}
END_TEST

START_TEST (test_insert_advance0)
{
    iproc_events_insert(events, 123);
    iproc_events_advance(events, 0);

    ck_assert_int_eq(iproc_events_ncur(events), 1);
    ck_assert_int_eq(iproc_events_npast(events), 0);
    ck_assert_int_eq(iproc_events_cur(events, 0), 123);
}
END_TEST

START_TEST (test_insert_advance)
{
    iproc_events_insert(events, 123);
    iproc_events_advance(events, 3);

    ck_assert_int_eq(iproc_events_ncur(events), 0);
    ck_assert_int_eq(iproc_events_npast(events), 1);
    ck_assert_int_eq(iproc_events_past(events, 0), 123);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 3);
}
END_TEST

START_TEST (test_insert_advance2)
{
    iproc_events_insert(events, 11);
    iproc_events_advance(events, 9);
    iproc_events_advance(events, 3);

    ck_assert_int_eq(iproc_events_ncur(events), 0);
    ck_assert_int_eq(iproc_events_npast(events), 1);
    ck_assert_int_eq(iproc_events_past(events, 0), 11);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 12);
}
END_TEST

START_TEST (test_insert3_advance)
{
    iproc_events_insert(events, 5);
    iproc_events_insert(events, 1);
    iproc_events_insert(events, 3);
    iproc_events_advance(events, 7);

    ck_assert_int_eq(iproc_events_ncur(events), 0);
    ck_assert_int_eq(iproc_events_npast(events), 3);
    ck_assert_int_eq(iproc_events_past(events, 0), 1);
    ck_assert_int_eq(iproc_events_past(events, 1), 3);
    ck_assert_int_eq(iproc_events_past(events, 2), 5);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 7);
    ck_assert_int_eq(iproc_events_past_dt(events, 1), 7);
    ck_assert_int_eq(iproc_events_past_dt(events, 2), 7);
}
END_TEST

START_TEST (test_insert_advance_insert_dup)
{
    iproc_events_insert(events, 1);
    iproc_events_insert(events, 3);
    iproc_events_advance(events, 10);
    iproc_events_insert(events, 1);
    iproc_events_advance(events, 5);

    ck_assert_int_eq(iproc_events_ncur(events), 0);
    ck_assert_int_eq(iproc_events_npast(events), 2);
    ck_assert_int_eq(iproc_events_past(events, 0), 1);
    ck_assert_int_eq(iproc_events_past(events, 1), 3);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 5);
    ck_assert_int_eq(iproc_events_past_dt(events, 1), 15);
}
END_TEST

START_TEST (test_insert_advance_complex)
{
    iproc_events_insert(events, 2);
    iproc_events_advance(events, 10);

    ck_assert_int_eq(iproc_events_npast(events), 1);
    ck_assert_int_eq(iproc_events_past(events, 0), 2);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 10);

    iproc_events_insert(events, 1);
    iproc_events_advance(events, 10);

    ck_assert_int_eq(iproc_events_npast(events), 2);
    ck_assert_int_eq(iproc_events_past(events, 0), 1);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 10);
    ck_assert_int_eq(iproc_events_past(events, 1), 2);
    ck_assert_int_eq(iproc_events_past_dt(events, 1), 20);

    iproc_events_insert(events, 4);
    iproc_events_advance(events, 10);

    ck_assert_int_eq(iproc_events_npast(events), 3);
    ck_assert_int_eq(iproc_events_past(events, 0), 1);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 20);
    ck_assert_int_eq(iproc_events_past(events, 1), 2);
    ck_assert_int_eq(iproc_events_past_dt(events, 1), 30);
    ck_assert_int_eq(iproc_events_past(events, 2), 4);
    ck_assert_int_eq(iproc_events_past_dt(events, 2), 10);

    iproc_events_insert(events, 3);
    iproc_events_advance(events, 10);

    ck_assert_int_eq(iproc_events_npast(events), 4);
    ck_assert_int_eq(iproc_events_past(events, 0), 1);
    ck_assert_int_eq(iproc_events_past_dt(events, 0), 30);
    ck_assert_int_eq(iproc_events_past(events, 1), 2);
    ck_assert_int_eq(iproc_events_past_dt(events, 1), 40);
    ck_assert_int_eq(iproc_events_past(events, 2), 3);
    ck_assert_int_eq(iproc_events_past_dt(events, 2), 10);
    ck_assert_int_eq(iproc_events_past(events, 3), 4);
    ck_assert_int_eq(iproc_events_past_dt(events, 3), 20);
}
END_TEST

START_TEST (test_find_past)
{
    int i;

    for (i = 0; i < 63; i++) {
        iproc_events_insert(events, i*2);
    }
    iproc_events_advance(events, 1);

    for (i = 0; i < 63; i++) {
        ck_assert_int_eq(iproc_events_find_past(events, i*2), i);
    }

    for (i = 0; i <= 63; i++) {
        ck_assert_int_eq(iproc_events_find_past(events, i*2 - 1), ~i);
    }
}
END_TEST

START_TEST (test_clear)
{
    iproc_events_insert(events, 2);
    iproc_events_advance(events, 7);
    iproc_events_insert(events, 1);

    iproc_events_clear(events);
    ck_assert_int_eq(iproc_events_ncur(events), 0);
    ck_assert_int_eq(iproc_events_npast(events), 0);
}
END_TEST

TCase *
events_tcase ()
{
    TCase *tc = tcase_create("events");
    tcase_add_checked_fixture(tc, setup, teardown);
    tcase_add_test(tc, test_init);
    tcase_add_test(tc, test_insert1);
    tcase_add_test(tc, test_insert2);
    tcase_add_test(tc, test_insert2_dup);
    tcase_add_test(tc, test_insert_advance0);
    tcase_add_test(tc, test_insert_advance);
    tcase_add_test(tc, test_insert_advance2);
    tcase_add_test(tc, test_insert3_advance);
    tcase_add_test(tc, test_insert_advance_insert_dup);
    tcase_add_test(tc, test_insert_advance_complex);
    tcase_add_test(tc, test_clear);
    tcase_add_test(tc, test_find_past);
    return tc;
}

Suite *
events_suite ()
{
    Suite *s = suite_create("events");
    suite_add_tcase(s, events_tcase());
    return s;
}
