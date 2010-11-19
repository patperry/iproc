#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string.h>
#include <glib.h>
#include <iproc/events.h>

static void
iproc_events_clear_cur (iproc_events *events)
{
    g_assert(events);
    events->ncur = 0;
}

static void
iproc_events_clear_past (iproc_events *events)
{
    g_assert(events);
    events->npast = 0;
}

static void
iproc_events_grow_cur (iproc_events *events)
{
    g_assert(events);
    g_assert(events->max_ncur > 0);

    int64_t max_ncur = 2 * events->max_ncur;
    events->cur = g_realloc(events->cur, sizeof(events->cur[0]) * max_ncur);
    events->max_ncur = max_ncur;
}

static void
iproc_events_grow_past (iproc_events *events)
{
    g_assert(events);
    g_assert(events->max_npast > 0);

    int64_t max_npast = 2 * events->max_npast;
    events->past    = g_realloc(events->past, sizeof(events->past[0]) * max_npast);
    events->past_dt = g_realloc(events->past_dt, sizeof(events->past_dt[0]) * max_npast);
    events->max_npast = max_npast;
}

static void
iproc_events_reserve_cur (iproc_events *events, int64_t n)
{
    g_assert(events);

    while (events->max_ncur < n) {
        iproc_events_grow_cur(events);
    }
}

static void
iproc_events_reserve_past (iproc_events *events, int64_t n)
{
    g_assert(events);

    while (events->max_npast < n) {
        iproc_events_grow_past(events);
    }
}

static void
iproc_events_move_past (iproc_events *events, int64_t start)
{
    g_assert(events);
    int64_t *tail    = events->past + start;
    int64_t *tail_dt = events->past_dt + start;
    int64_t n = events->npast;
    int64_t ntail = n - start;

    g_memmove(tail + 1, tail, ntail * sizeof(tail[0]));
    g_memmove(tail_dt + 1, tail_dt, ntail * sizeof(tail_dt[0]));
}

iproc_events *
iproc_events_new ()
{
    iproc_events *events = g_malloc(sizeof(*events));

    g_return_val_if_fail(events, NULL);

    events->ncur = 0;
    events->max_ncur = 1;
    events->cur = g_malloc(sizeof(events->cur[0]) * events->max_ncur);
    events->npast = 0;
    events->max_npast = 1;
    events->past = g_malloc(sizeof(events->past[0]) * events->max_npast);
    events->past_dt = g_malloc(sizeof(events->past_dt[0]) * events->max_npast);

    if (!(events->cur && events->past && events->past_dt)) {
        iproc_events_free(events);
        events = NULL;
    }

    return events;
}

void
iproc_events_free (iproc_events *events)
{
    if (events) {
        g_free(events->cur);
        g_free(events->past);
        g_free(events->past_dt);
        g_free(events);
    }
}

void
iproc_events_clear (iproc_events *events)
{
    iproc_events_clear_cur(events);
    iproc_events_clear_past(events);
}

void
iproc_events_insert (iproc_events *events,
                     int64_t       e)
{
    g_assert(events);

    if (iproc_events_find_cur(events, e) < 0) {
        int64_t ncur = events->ncur;

        iproc_events_reserve_cur(events, ncur + 1);

        events->cur[ncur] = e;
        events->ncur = ncur + 1;
    }

    g_assert(events->ncur <= events->max_ncur);
}

void
iproc_events_advance (iproc_events *events,
                      int64_t       dt)
{
    g_assert(events);
    g_assert(dt >= 0);

    int64_t e;
    int64_t ic, ip;
    int64_t nc = iproc_events_ncur(events);
    int64_t np = iproc_events_npast(events);

    if (dt <= 0)
        return;

    /* Add time to all past events */
    for (ip = 0; ip < np; ip++) {
        events->past_dt[ip] = events->past_dt[ip] + dt;
    }

    /* Move current events to past event set */
    for (ic = 0; ic < nc; ic++) {
        e = iproc_events_cur(events, ic);
        ip = iproc_events_find_past(events, e);

        /* If event doesn't already exist in past set, insert it */
        if (ip < 0) {
            np++;
            ip = ~ip;
            iproc_events_reserve_past(events, np);
            iproc_events_move_past(events, ip);
            events->past[ip] = e;
            events->npast = np;
        }

        /* Set the time of the event to the time advance */
        events->past_dt[ip] = dt;
    }

    iproc_events_clear_cur(events);
}

int64_t
iproc_events_ncur (iproc_events *events)
{
    g_assert(events);
    return events->ncur;
}

int64_t
iproc_events_cur (iproc_events *events,
                  int64_t       i)
{
    g_assert(events);
    g_assert(0 <= i);
    g_assert(i < events->ncur);

    return events->cur[i];
}

int64_t
iproc_events_find_cur (iproc_events *events,
                       int64_t e)
{
    g_assert(events);
    int64_t i, n = events->ncur;

    for (i = 0; i < n; i++) {
        if (iproc_events_cur(events, i) == e)
            break;
    }

    if (i == n)
        i = ~i;

    return i;
}

static int64_t
iproc_lsearch (int64_t n, int64_t *array, int64_t key)
{
    g_assert(n >= 0);
    g_assert(array);

    int64_t i;
    int64_t x;

    for (i = 0; i < n; i++) {
        x = array[i];

        if (key == x) {
            break;
        } else if (key < x) {
            i = ~i;
            break;
        }
    }

    if (i == n)
        i = ~i;

    return i;
}

static int64_t
iproc_bsearch (int64_t n, int64_t *array, int64_t key)
{
    g_assert(n >= 0);
    g_assert(array);

    int64_t begin = 0;
    int64_t end = n;
    int64_t target;
    int64_t x;

    while (begin < end) {
        target = begin + ((end - begin) >> 1);
        x = array[target];
        if (key < x) {
            end = target;
        } else if (key > x) {
            begin = target + 1;
        } else { /* key == x */
            return target;
        }
    }

    return ~end; /* begin == end, not found */
}

int64_t
iproc_events_find_past (iproc_events *events, int64_t e)
{
    g_assert(events);
    return iproc_bsearch(events->npast, events->past, e);
}

int64_t
iproc_events_npast (iproc_events *events)
{
    g_assert(events);
    return events->npast;
}

int64_t
iproc_events_past (iproc_events *events,
                   int64_t        i)
{
    g_assert(events);
    g_assert(0 <= i);
    g_assert(i < iproc_events_npast(events));

    return events->past[i];
}

int64_t
iproc_events_past_dt (iproc_events *events,
                      int64_t        i)
{
    g_assert(events);
    g_assert(0 <= i);
    g_assert(i < iproc_events_npast(events));

    return events->past_dt[i];
}
