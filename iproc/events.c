#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <string.h>
#include <assert.h>
#include <iproc/events.h>
#include <iproc/memory.h>

static int
compare_event (void *px, void *py)
{
    int64_t x = ((iproc_event *) px)->e;
    int64_t y = ((iproc_event *) py)->e;

    if (x < y) {
        return -1;
    } else if (x > y) {
        return 1;
    } else {
        return 0;
    }
}

static void
iproc_events_clear_cur (iproc_events *events)
{
    assert(events);
    iproc_array_set_size(events->cur, 0);
}

static void
iproc_events_clear_past (iproc_events *events)
{
    assert(events);
    iproc_array_set_size(events->past, 0);
}

iproc_events *
iproc_events_new ()
{
    iproc_events *events = iproc_malloc(sizeof(*events));

    if (!events) return NULL;

    events->cur  = iproc_array_new(sizeof(iproc_event));
    events->past = iproc_array_new(sizeof(iproc_past_event));

    if (!(events->cur && events->past)) {
        iproc_events_free(events);
        events = NULL;
    }

    return events;
}

void
iproc_events_free (iproc_events *events)
{
    if (events) {
        iproc_array_unref(events->cur);
        iproc_array_unref(events->past);
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
    iproc_event event = { e };

    if (iproc_events_find_cur(events, e) < 0) {
        iproc_array_append(events->cur, &event);
    }
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
    iproc_array *past = events->past;

    if (dt <= 0)
        return;

    /* Add time to all past events */
    for (ip = 0; ip < np; ip++) {
        iproc_array_index(past, iproc_past_event, ip).dt += dt;
    }

    /* Move current events to past event set */
    for (ic = 0; ic < nc; ic++) {
        e = iproc_events_cur(events, ic);
        ip = iproc_events_find_past(events, e);

        /* If event doesn't already exist in past set, insert it */
        if (ip < 0) {
            iproc_event event = { e=e };
            iproc_past_event past_event = { event=event, dt=dt };
            iproc_array_insert(past, ~ip, &past_event);

        /* Set the time of the event to the time advance */
        } else {
            iproc_array_index(past, iproc_past_event, ip).dt = dt;
        }
    }

    iproc_events_clear_cur(events);
}

int64_t
iproc_events_ncur (iproc_events *events)
{
    assert(events);
    return iproc_array_size(events->cur);
}

int64_t
iproc_events_cur (iproc_events *events,
                  int64_t       i)
{
    assert(events);
    assert(0 <= i);
    assert(i < iproc_events_ncur(events));

    return iproc_array_index(events->cur, iproc_event, i).e;
}

int64_t
iproc_events_find_cur (iproc_events *events,
                       int64_t       e)
{
    assert(events);
    iproc_event event = { e=e };
    return iproc_array_lfind(events->cur, &event, compare_event);
}

int64_t
iproc_events_find_past (iproc_events *events,
                        int64_t       e)
{
    assert(events);

    /* use that iproc_past_event inherits from iproc_event */
    iproc_event event = { e=e };
    return iproc_array_bsearch(events->past, &event, compare_event);
}

int64_t
iproc_events_npast (iproc_events *events)
{
    assert(events);
    return iproc_array_size(events->past);
}

int64_t
iproc_events_past (iproc_events *events,
                   int64_t       i)
{
    assert(events);
    assert(0 <= i);
    assert(i < iproc_events_npast(events));

    return iproc_array_index(events->past, iproc_past_event, i).event.e;
}

int64_t
iproc_events_past_dt (iproc_events *events,
                      int64_t       i)
{
    assert(events);
    assert(0 <= i);
    assert(i < iproc_events_npast(events));

    return iproc_array_index(events->past, iproc_past_event, i).dt;
}
