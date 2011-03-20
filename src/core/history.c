#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "history.h"

static void
trace_array_grow (iproc_array *array,
                  int64_t      n)
{
    assert(array);
    int64_t nold = iproc_array_size(array);

    if (n > nold) {
        iproc_array_set_size(array, n);
    }
}

static void
trace_array_clear (iproc_array *array)
{
    assert(array);
    int64_t n = iproc_array_size(array);
    int64_t i;
    iproc_history_trace *ht;

    for (i = 0; i < n; i++) {
        ht = &(iproc_array_index(array, iproc_history_trace, i));
        ht->tcur = -INFINITY;

        if (ht->trace)
            iproc_trace_clear(ht->trace);
    }
}

static iproc_trace *
trace_array_get (double       tcur,
                 iproc_array *array,
                 int64_t      i)
{
    assert(array);
    assert(i >= 0);

    trace_array_grow(array, i + 1);

    iproc_history_trace *ht = &iproc_array_index(array,
                                                 iproc_history_trace,
                                                 i);
    iproc_trace *t;

    if (!(ht->trace)) {
        ht->trace = iproc_trace_new(tcur);
    }

    t = ht->trace;

    if (ht->tcur != tcur) {
        iproc_trace_advance_to(t, tcur);
        ht->tcur = tcur;
    }

    return t;
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
    history->send = iproc_array_new(sizeof(iproc_history_trace));
    history->recv = iproc_array_new(sizeof(iproc_history_trace));
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
    trace_array_clear(history->send);
    trace_array_clear(history->recv);
}

double
iproc_history_tcur (iproc_history *history)
{
    assert(history);
    return history->tcur;
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

    iproc_trace *efrom = iproc_history_send(history, from);
    iproc_trace *eto   = iproc_history_recv(history, to);

    iproc_trace_insert(efrom, to);
    iproc_trace_insert(eto, from);
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

iproc_trace *
iproc_history_send (iproc_history *history,
                    int64_t        i)
{
    assert(history);
    assert(0 <= i);

    return trace_array_get(history->tcur, history->send, i);
}

iproc_trace *
iproc_history_recv (iproc_history *history,
                    int64_t        j)
{
    assert(history);
    assert(0 <= j);

    return trace_array_get(history->tcur, history->recv, j);
}
