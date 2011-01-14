#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include "memory.h"
#include "messages.h"

static void
iproc_messages_free (iproc_messages *msgs)
{
    if (msgs) {
        iproc_array_unref(msgs->array);
        iproc_array_unref(msgs->recipients);
        iproc_free(msgs);
    }
}

iproc_messages *
iproc_messages_new (int64_t t0)
{
    iproc_messages *msgs = iproc_malloc(sizeof(*msgs));
    if (!msgs)
        return NULL;

    msgs->tcur = t0;
    msgs->array = iproc_array_new(sizeof(iproc_message));
    msgs->recipients = iproc_array_new(sizeof(int64_t));
    msgs->max_to = -1;
    msgs->max_from = -1;
    msgs->max_nto = 0;
    iproc_refcount_init(&msgs->refcount);

    if (!(msgs->array && msgs->recipients)) {
        iproc_messages_free(msgs);
        msgs = NULL;
    }

    return msgs;
}

iproc_messages *
iproc_messages_ref (iproc_messages *msgs)
{
    if (msgs) {
        iproc_refcount_get(&msgs->refcount);
    }
    return msgs;
}

static void
iproc_messages_release (iproc_refcount *refcount)
{
    iproc_messages *msgs = container_of(refcount, iproc_messages, refcount);
    iproc_messages_free(msgs);
}

iproc_messages *
iproc_messages_unref (iproc_messages *msgs)
{
    if (msgs) {
        iproc_refcount_put(&msgs->refcount, iproc_messages_release);
    }
}

int64_t
iproc_messages_size (iproc_messages *msgs)
{
    assert(msgs);
    return iproc_array_size(msgs->array);
}

void
iproc_messages_advance (iproc_messages *msgs,
                        int64_t         dt)
{
    assert(msgs);
    assert(dt >= 0);
    msgs->tcur += dt;
}

void
iproc_messages_advance_to (iproc_messages *msgs,
                           int64_t         t)
{
    assert(msgs);
    assert(t >= msgs->tcur);
    msgs->tcur = t;
}

void
iproc_messages_insert (iproc_messages *msgs,
                       int64_t         from,
                       int64_t         to)
{
    assert(msgs);
    assert(from >= 0);
    assert(to >= 0);
    iproc_messages_insertm(msgs, from, &to, 1);
}

void
iproc_messages_insertm (iproc_messages *msgs,
                        int64_t         from,
                        int64_t        *to,
                        int64_t         nto)
{
    assert(msgs);
    assert(from >= 0);
    assert(nto >= 0);
    assert(to || nto == 0);

    int64_t time = msgs->tcur;
    iproc_array *array = msgs->array;
    iproc_array *recipients = msgs->recipients;

    int64_t n = iproc_array_size(recipients);
    int64_t *mto  = &(iproc_array_index(recipients, int64_t, n));
    iproc_message m = {
        time = time, 
        from = from,
        to  = mto,
        nto = nto
    };

    int64_t i;
    for (i = 0; i < nto; i++) {
        assert(to[i] >= 0);

        iproc_array_append(recipients, to + i);
        if (to[i] > msgs->max_to)
            msgs->max_to = to[i];
    }

    if (from > msgs->max_from)
        msgs->max_from = from;
    if (nto > msgs->max_nto)
        msgs->max_nto = nto;

    iproc_array_append(array, &m);
}

int64_t
iproc_messages_max_from (iproc_messages *msgs)
{
    assert(msgs);
    return msgs->max_from;
}

int64_t
iproc_messages_max_to (iproc_messages *msgs)
{
    assert(msgs);
    return msgs->max_to;
}

int64_t
iproc_messages_max_nto (iproc_messages *msgs)
{
    assert(msgs);
    return msgs->max_nto;
}

