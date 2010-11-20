#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string.h>
#include <assert.h>
#include <iproc/events.h>
#include <iproc/memory.h>

static void
iproc_events_clear_cur (iproc_events *events)
{
    assert(events);
    events->ncur = 0;
}

static void
iproc_events_clear_past (iproc_events *events)
{
    assert(events);
    events->npast = 0;
}

static void
iproc_events_grow_cur (iproc_events *events)
{
    assert(events);
    assert(events->max_ncur > 0);

    int64_t max_ncur = 2 * events->max_ncur;
    events->cur = iproc_realloc(events->cur, sizeof(events->cur[0]) * max_ncur);
    events->max_ncur = max_ncur;
}

static void
iproc_events_grow_past (iproc_events *events)
{
    assert(events);
    assert(events->max_npast > 0);

    int64_t max_npast = 2 * events->max_npast;
    events->past    = iproc_realloc(events->past, sizeof(events->past[0]) * max_npast);
    events->past_dt = iproc_realloc(events->past_dt, sizeof(events->past_dt[0]) * max_npast);
    events->max_npast = max_npast;
}

static void
iproc_events_reserve_cur (iproc_events *events,
                          int64_t       n)
{
    assert(events);

    while (events->max_ncur < n) {
        iproc_events_grow_cur(events);
    }
}

static void
iproc_events_reserve_past (iproc_events *events,
                           int64_t       n)
{
    assert(events);

    while (events->max_npast < n) {
        iproc_events_grow_past(events);
    }
}

static void
iproc_events_move_past (iproc_events *events,
                        int64_t       start)
{
    assert(events);
    int64_t *tail    = events->past + start;
    int64_t *tail_dt = events->past_dt + start;
    int64_t n = events->npast;
    int64_t ntail = n - start;

    memmove(tail + 1, tail, ntail * sizeof(tail[0]));
    memmove(tail_dt + 1, tail_dt, ntail * sizeof(tail_dt[0]));
}

iproc_events *
iproc_events_new ()
{
    iproc_events *events = iproc_malloc(sizeof(*events));

    if (!events) return NULL;

    events->ncur = 0;
    events->max_ncur = 1;
    events->cur = iproc_malloc(sizeof(events->cur[0]) * events->max_ncur);
    events->npast = 0;
    events->max_npast = 1;
    events->past = iproc_malloc(sizeof(events->past[0]) * events->max_npast);
    events->past_dt = iproc_malloc(sizeof(events->past_dt[0]) * events->max_npast);

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
        iproc_free(events->cur);
        iproc_free(events->past);
        iproc_free(events->past_dt);
        iproc_free(events);
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
    assert(events);

    if (iproc_events_find_cur(events, e) < 0) {
        int64_t ncur = events->ncur;

        iproc_events_reserve_cur(events, ncur + 1);

        events->cur[ncur] = e;
        events->ncur = ncur + 1;
    }

    assert(events->ncur <= events->max_ncur);
}

void
iproc_events_advance (iproc_events *events,
                      int64_t       dt)
{
    assert(events);
    assert(dt >= 0);

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
    assert(events);
    return events->ncur;
}

int64_t
iproc_events_cur (iproc_events *events,
                  int64_t       i)
{
    assert(events);
    assert(0 <= i);
    assert(i < events->ncur);

    return events->cur[i];
}

int64_t
iproc_events_find_cur (iproc_events *events,
                       int64_t       e)
{
    assert(events);
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
iproc_lsearch (int64_t n,
               int64_t *array,
               int64_t key)
{
    assert(n >= 0);
    assert(array);

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
iproc_bsearch (int64_t n,
               int64_t *array,
               int64_t key)
{
    assert(n >= 0);
    assert(array);

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
iproc_events_find_past (iproc_events *events,
                        int64_t       e)
{
    assert(events);
    int64_t  n = events->npast;
    int64_t *array = events->past;
    int64_t  key = e;
    int64_t  i;
    
    if (n <= 8) {
        i = iproc_lsearch(n, array, key);
    } else {
        i = iproc_bsearch(n, array, key);
    }
    
    return i;
}

int64_t
iproc_events_npast (iproc_events *events)
{
    assert(events);
    return events->npast;
}

int64_t
iproc_events_past (iproc_events *events,
                   int64_t       i)
{
    assert(events);
    assert(0 <= i);
    assert(i < iproc_events_npast(events));

    return events->past[i];
}

int64_t
iproc_events_past_dt (iproc_events *events,
                      int64_t       i)
{
    assert(events);
    assert(0 <= i);
    assert(i < iproc_events_npast(events));

    return events->past_dt[i];
}
