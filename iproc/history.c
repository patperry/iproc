#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <iproc/memory.h>
#include <iproc/history.h>

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
        he->elapsed = 0;

        if (he->events)
            iproc_events_clear(he->events);
    }
}

static iproc_events *
iproc_history_get (int64_t      elapsed,
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
        he->events = iproc_events_new();
    }

    e = he->events;

    if (he->elapsed != elapsed) {
        iproc_events_advance(e, elapsed - he->elapsed);
        he->elapsed = elapsed;
    }

    return e;
}

iproc_history *
iproc_history_new ()
{
    iproc_history *history = iproc_malloc(sizeof(*history));

    if (!history) return NULL;

    history->send = iproc_array_new(sizeof(iproc_history_events));
    history->recv = iproc_array_new(sizeof(iproc_history_events));

    if (!(history->send && history->recv)) {
        iproc_history_free(history);
        history = NULL;
    }

    return history;
}

void
iproc_history_free (iproc_history *history)
{
    if (history) {
        iproc_array_free(history->send);
        iproc_array_free(history->recv);
        iproc_free(history);
    }
}

void
iproc_history_clear (iproc_history *history)
{
    assert(history);
    history->elapsed = 0;
    iproc_history_clear_events(history->send);
    iproc_history_clear_events(history->recv);
}

void
iproc_history_advance (iproc_history *history,
                       int64_t        dt)
{
    assert(history);
    assert(dt >= 0);

    if (dt > 0) {
        history->elapsed += dt;
    }
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
                       int            nfrom,
                       int64_t       *from,
                       int            nto,
                       int64_t       *to)
{
    assert(history);
    assert(nfrom > 0);
    assert(from);
    assert(nto > 0);
    assert(to);

    int i, j;

    for (i = 0; i < nfrom; i++) {
        for (j = 0; j < nto; j++) {
            iproc_history_insert(history, from[i], to[j]);
        }
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

    return iproc_history_get(history->elapsed, history->send, i);
}

iproc_events *
iproc_history_recv (iproc_history *history,
                    int64_t        j)
{
    assert(history);
    assert(0 <= j);

    return iproc_history_get(history->elapsed, history->recv, j);
}
