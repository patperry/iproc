#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "history.h"

static void
iproc_history_grow_events (iproc_array *array,
                           int64_t      n)
{
    assert(array);
    int64_t nold = iproc_array_size(array);

    if (n > nold) {
        iproc_array_set_size(array, n);
    }
}

static void
iproc_history_clear_events (iproc_array *array)
{
    assert(array);
    int64_t n = iproc_array_size(array);
    int64_t i;
    iproc_history_events *he;

    for (i = 0; i < n; i++) {
        he = &(iproc_array_index(array, iproc_history_events, i));
        he->tcur = -INFINITY;

        if (he->events)
            iproc_events_clear(he->events);
    }
}

static iproc_events *
iproc_history_get (double       tcur,
                   iproc_array *array,
                   int64_t      i)
{
    assert(array);
    assert(i >= 0);

    iproc_history_grow_events(array, i + 1);

    iproc_history_events *he = &(iproc_array_index(array,
                                                   iproc_history_events,
                                                   i));
    iproc_events *e;

    if (!(he->events)) {
        he->events = iproc_events_new(tcur);
    }

    e = he->events;

    if (he->tcur != tcur) {
        iproc_events_advance_to(e, tcur);
        he->tcur = tcur;
    }

    return e;
}

static void
iproc_history_free (iproc_history *history)
{
    if (history) {
        iproc_array_unref(history->send);
        iproc_array_unref(history->recv);
        iproc_free(history);
    }
}

iproc_history *
iproc_history_new ()
{
    iproc_history *history = iproc_malloc(sizeof(*history));

    if (!history) return NULL;

    history->tcur = -INFINITY;
    history->send = iproc_array_new(sizeof(iproc_history_events));
    history->recv = iproc_array_new(sizeof(iproc_history_events));
    iproc_refcount_init(&history->refcount);

    if (!(history->send && history->recv)) {
        iproc_history_free(history);
        history = NULL;
    }

    return history;
}

iproc_history *
iproc_history_ref (iproc_history *history)
{
    if (history) {
        iproc_refcount_get(&history->refcount);
    }
    return history;
}

static void
iproc_history_release (iproc_refcount *refcount)
{
    iproc_history *history = container_of(refcount, iproc_history, refcount);
    iproc_history_free(history);
}

void
iproc_history_unref (iproc_history *history)
{
    if (!history)
        return;

    iproc_refcount_put(&history->refcount, iproc_history_release);
}

void
iproc_history_clear (iproc_history *history)
{
    assert(history);
    history->tcur = -INFINITY;
    iproc_history_clear_events(history->send);
    iproc_history_clear_events(history->recv);
}

void
iproc_history_advance_to (iproc_history *history,
                          double         t)
{
    assert(history);
    assert(history->tcur <= t);

    history->tcur = t;
}

void
iproc_history_insert (iproc_history *history,
                      int64_t        from,
                      int64_t        to)
{
    assert(history);
    assert(from >= 0);
    assert(to >= 0);

    iproc_events *efrom = iproc_history_send(history, from);
    iproc_events *eto   = iproc_history_recv(history, to);

    iproc_events_insert(efrom, to);
    iproc_events_insert(eto, from);
}

void
iproc_history_insertm (iproc_history *history,
                       int64_t        from,
                       int64_t       *to,
                       int64_t        nto)
{
    assert(history);
    assert(to || nto == 0);
    assert(nto >= 0);

    int i;

    for (i = 0; i < nto; i++) {
        iproc_history_insert(history, from, to[i]);
    }
}

int64_t
iproc_history_nsend (iproc_history *history)
{
    assert(history);
    return iproc_array_size(history->send);
}

int64_t
iproc_history_nrecv (iproc_history *history)
{
    assert(history);
    return iproc_array_size(history->recv);
}

iproc_events *
iproc_history_send (iproc_history *history,
                    int64_t        i)
{
    assert(history);
    assert(0 <= i);

    return iproc_history_get(history->tcur, history->send, i);
}

iproc_events *
iproc_history_recv (iproc_history *history,
                    int64_t        j)
{
    assert(history);
    assert(0 <= j);

    return iproc_history_get(history->tcur, history->recv, j);
}