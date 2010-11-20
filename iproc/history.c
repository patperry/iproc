#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <iproc/memory.h>
#include <iproc/history.h>


static void
iproc_history_reserve_events (int64_t               *pn,
                              iproc_history_events **parray,
                              int64_t                n)
{
    assert(pn);
    assert(parray);

    iproc_history_events *array = *parray;
    int64_t nold = *pn;
    int64_t i;

    if (n > nold) {
        array = iproc_realloc(array, n * sizeof(*array));
        for (i = nold; i < n; i++) {
            array[i].elapsed = 0;
            array[i].events = iproc_events_new();
        }

        *pn = n;
        *parray = array;
    }
}

static void
iproc_history_grow_events (int64_t               *pn,
                           int64_t               *pmaxn,
                           iproc_history_events **parray,
                           int64_t                n)
{
    assert(pn);
    assert(pmaxn);
    assert(parray);

    if (n > *pn) {
        iproc_history_reserve_events(pmaxn, parray, n);
        *pn = n;
    }
}

static void
iproc_history_free_events (int64_t n, iproc_history_events *array)
{
    assert (n >= 0);
    int i;

    if (array) {
        for (i = 0; i < n; i++) {
            iproc_events_free(array[i].events);
        }
        iproc_free(array);
    }
}

static void
iproc_history_clear_events (int64_t n, iproc_history_events *array)
{
    assert (n >= 0);
    int i;

    if (array) {
        for (i = 0; i < n; i++) {
            array[i].elapsed = 0;
            iproc_events_clear(array[i].events);
        }
    }
}

static iproc_events *
iproc_history_get (int64_t               elapsed,
                   iproc_history_events *array,
                   int64_t               i)
{
    assert(array);
    assert(i >= 0);

    iproc_history_events *he = array + i;
    iproc_events *e = he->events;
    int64_t he_elapsed = he->elapsed;

    if (he_elapsed != elapsed) {
        iproc_events_advance(e, elapsed - he_elapsed);
        he->elapsed = elapsed;
    }

    return e;
}

static void
iproc_history_grow_send (iproc_history *history,
                         int64_t        nsend)
{
    assert(history);
    assert(nsend >= 0);

    iproc_history_grow_events(&(history->nsend),
                              &(history->max_nsend),
                              &(history->send),
                              nsend);
}

static void
iproc_history_grow_recv (iproc_history *history,
                         int64_t        nrecv)
{
    assert(history);
    assert(nrecv >= 0);

    iproc_history_grow_events(&(history->nrecv),
                              &(history->max_nrecv),
                              &(history->recv),
                              nrecv);
}

iproc_history *
iproc_history_new ()
{
    iproc_history *history = iproc_malloc(sizeof(*history));

    if (!history) return NULL;

    history->send = NULL;
    history->recv = NULL;
    history->elapsed = 0;
    history->nsend = 0;
    history->nrecv = 0;
    history->max_nsend = 0;
    history->max_nrecv = 0;
    return history;
}

void
iproc_history_free (iproc_history *history)
{
    if (history) {
        iproc_history_free_events(history->max_nsend, history->send);
        iproc_history_free_events(history->max_nrecv, history->recv);
        iproc_free(history);
    }
}

void
iproc_history_clear (iproc_history *history)
{
    assert(history);
    history->elapsed = 0;
    history->nsend = 0;
    history->nrecv = 0;
    iproc_history_clear_events(history->max_nsend, history->send);
    iproc_history_clear_events(history->max_nrecv, history->recv);
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

    iproc_history_grow_send(history, from + 1);
    iproc_history_grow_recv(history, to + 1);

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
    return history->nsend;
}

int64_t
iproc_history_nrecv (iproc_history *history)
{
    assert(history);
    return history->nrecv;
}

iproc_events *
iproc_history_send (iproc_history *history,
                    int64_t        i)
{
    assert(history);
    assert(0 <= i);
    assert(i < iproc_history_nsend(history));

    return iproc_history_get(history->elapsed, history->send, i);
}

iproc_events *
iproc_history_recv (iproc_history *history,
                    int64_t        j)
{
    assert(history);
    assert(0 <= j);
    assert(j < iproc_history_nrecv(history));

    return iproc_history_get(history->elapsed, history->recv, j);
}
