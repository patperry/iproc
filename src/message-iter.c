
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

int64_t
iproc_message_iter_time (iproc_message_iter *it)
{
    assert(it);
    assert(it->message);
    return it->message->time;
}

int64_t
iproc_message_iter_nties (iproc_message_iter *it)
{
    assert(it);
    assert(it->message);
    return it->nties;
}

void
iproc_message_iter_select (iproc_message_iter *it,
                           int64_t             tie)
{
    assert(it);
    assert(tie >= 0);
    assert(tie < iproc_message_iter_nties(it));

    int64_t i = it->offset + tie;
    it->message = &(iproc_array_index(it->messages->array, iproc_message, i));
}

int64_t
iproc_message_iter_from (iproc_message_iter *it)
{
    assert(it);
    assert(it->message);
    return it->message->from;
}

int64_t
iproc_message_iter_nto (iproc_message_iter *it)
{
    assert(it);
    assert(it->message);
    return it->message->nto;
}

int64_t *
iproc_message_iter_to (iproc_message_iter *it)
{
    assert(it);
    assert(it->message);
    return it->message->to;
}

void
iproc_message_iter_reset (iproc_message_iter *it)
{
    assert(it);
    it->message = NULL;
    it->offset = 0;
    it->nties = 0;
}

int
iproc_message_iter_next (iproc_message_iter *it)
{
    assert(it);
    int64_t offset = it->offset + it->nties;
    int64_t nties = 0;
    iproc_message *message = NULL;

    iproc_array *messages = it->messages->array;
    int64_t n = iproc_array_size(messages);
    int64_t has_next = offset < n ? 1 : 0;

    if (has_next) {
        message = &(iproc_array_index(messages, iproc_message, offset));
        int64_t time = message[0].time;
        int64_t nties_max = n - offset;
        nties = 1;

        while (nties < nties_max && message[nties].time == time) {
            nties++;
        }
    }

    it->offset = offset;
    it->nties = nties;
    it->message = message;
    return has_next;
}
