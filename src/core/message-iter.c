
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include "memory.h"
#include "messages.h"

static void
iproc_message_iter_free (iproc_message_iter *it)
{
    if (it) {
        iproc_history_unref(it->history);
        iproc_messages_unref(it->messages);
        iproc_free(it);
    }
}

iproc_message_iter *
iproc_message_iter_new (iproc_messages *msgs)
{
    iproc_message_iter *it = iproc_malloc(sizeof(*it));
    if (!it)
        return NULL;

    it->messages = iproc_messages_ref(msgs);
    it->history = iproc_history_new();
    iproc_refcount_init(&it->refcount);
    iproc_message_iter_reset(it);

    return it;
}

iproc_message_iter *
iproc_message_iter_ref (iproc_message_iter *it)
{
    if (it) {
        iproc_refcount_get(&it->refcount);
    }
    
    return it;
}

static void
iproc_message_iter_release (iproc_refcount *refcount)
{
    iproc_message_iter *it = container_of(refcount, iproc_message_iter, refcount);
    iproc_message_iter_free(it);
}

void
iproc_message_iter_unref (iproc_message_iter *it)
{
    if (it) {
        iproc_refcount_put(&it->refcount, iproc_message_iter_release);
    }
}

double
iproc_message_iter_time (iproc_message_iter *it)
{
    assert(iproc_message_iter_started(it));
    assert(!iproc_message_iter_finished(it));
    return it->message->time;
}

int64_t
iproc_message_iter_ntie (iproc_message_iter *it)
{
    assert(iproc_message_iter_started(it));
    assert(!iproc_message_iter_finished(it));
    return it->ntie;
}

void
iproc_message_iter_select (iproc_message_iter *it,
                           int64_t             tie)
{
    assert(iproc_message_iter_started(it));
    assert(!iproc_message_iter_finished(it));
    assert(tie >= 0);
    assert(tie < iproc_message_iter_ntie(it));

    int64_t i = it->offset + tie;
    it->message = &(iproc_array_index(it->messages->array, iproc_message, i));
}

int64_t
iproc_message_iter_from (iproc_message_iter *it)
{
    assert(iproc_message_iter_started(it));
    assert(!iproc_message_iter_finished(it));

    return it->message->from;
}

int64_t
iproc_message_iter_nto (iproc_message_iter *it)
{
    assert(iproc_message_iter_started(it));
    assert(!iproc_message_iter_finished(it));

    return it->message->nto;
}

int64_t *
iproc_message_iter_to (iproc_message_iter *it)
{
    assert(iproc_message_iter_started(it));
    assert(!iproc_message_iter_finished(it));

    int64_t ito = it->message->ito;
    int64_t *to = &(iproc_array_index(it->messages->recipients, int64_t, ito));
    return to;
}

void
iproc_message_iter_reset (iproc_message_iter *it)
{
    if (!it)
        return;

    iproc_history_clear(it->history);
    it->message = NULL;
    it->offset = 0;
    it->ntie = 0;
    it->finished = 0;
}

int
iproc_message_iter_next (iproc_message_iter *it)
{
    if (iproc_message_iter_finished(it))
        return 0;

    int64_t offset = it->offset + it->ntie;

    iproc_array *messages = it->messages->array;
    iproc_array *recipients = it->messages->recipients;
    iproc_history *history = it->history;
    int64_t n = iproc_array_size(messages);
    int has_next = offset < n ? 1 : 0;

    if (has_next) {
        iproc_message *message = &(iproc_array_index(messages, iproc_message, offset));
        double time = message[0].time;
        int64_t ntie_max = n - offset;
        int64_t ntie = 0;

        iproc_history_advance_to(history, time);

        while (ntie < ntie_max && message[ntie].time == time) {
            int64_t  msg_from = message[ntie].from;
            int64_t  msg_ito = message[ntie].ito;
            int64_t *msg_to = &(iproc_array_index(recipients, int64_t, msg_ito));
            int64_t  msg_nto = message[ntie].nto;

            iproc_history_insertm(history, msg_from, msg_to, msg_nto);
            ntie++;
        }
        it->ntie = ntie;
        it->message = message;
    } else {
        it->finished = 1;
    }
    it->offset = offset;

    return has_next;
}

int
iproc_message_iter_started (iproc_message_iter *it)
{
    if (!it)
        return 0;

    return it->message != NULL;
}

int
iproc_message_iter_finished (iproc_message_iter *it)
{
    if (!it)
        return 1;

    return it->finished;
}

iproc_history *
iproc_message_iter_history (iproc_message_iter *it)
{
    if (!it)
        return NULL;

    return it->history;
}
