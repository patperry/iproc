#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include "memory.h"
#include "cursor.h"

static void
iproc_cursor_free (iproc_cursor *cursor)
{
    if (cursor) {
        iproc_message_iter_unref(cursor->it);
        iproc_history_unref(cursor->history);
        iproc_free(cursor);
    }
}

iproc_cursor *
iproc_cursor_new (iproc_messages *msgs)
{
    iproc_cursor *cursor = iproc_malloc(sizeof(*cursor));

    if (!cursor)
        return NULL;

    cursor->it = iproc_message_iter_new(msgs);
    cursor->history = iproc_history_new();
    iproc_refcount_init(&cursor->refcount);

    if (!(cursor->it && cursor->history)) {
        iproc_cursor_free(cursor);
        cursor = NULL;
    } else {
        iproc_cursor_reset(cursor);
    }

    return cursor;
}

iproc_cursor *
iproc_cursor_ref (iproc_cursor *cursor)
{
    if (cursor) {
        iproc_refcount_get(&cursor->refcount);
    }
    return cursor;
}

static void
iproc_cursor_release (iproc_refcount *refcount)
{
    iproc_cursor *cursor = container_of(refcount, iproc_cursor, refcount);
    iproc_cursor_free(cursor);
}

void
iproc_cursor_unref (iproc_cursor *cursor)
{
    if (cursor) {
        iproc_refcount_put(&cursor->refcount, iproc_cursor_release);
    }
}
    
void
iproc_cursor_reset (iproc_cursor *cursor)
{
    if (!cursor)
        return;

    cursor->tcur = INT64_MIN;
    iproc_history_clear(cursor->history);
    iproc_message_iter_reset(cursor->it);
}

int
iproc_cursor_next (iproc_cursor *cursor)
{
    if (!cursor)
        return 0;

    int64_t tprev = iproc_cursor_time(cursor);
    iproc_message_iter *it = cursor->it;
    int has_next = 0;

    if (iproc_message_iter_next(it)) {
        iproc_history *history = cursor->history;
        int64_t tnext = iproc_message_iter_time(it);
        int64_t delta = tnext - tprev;
        int64_t i, n = iproc_message_iter_ntie(it);

        iproc_history_advance(history, delta);
        for (i = 0; i < n; i++) {
            iproc_message_iter_select(it, i);

            int64_t  msg_from = iproc_message_iter_from(it);
            int64_t *msg_to   = iproc_message_iter_to(it);
            int64_t  msg_nto  = iproc_message_iter_nto(it);

            iproc_history_insertm(history, msg_from, msg_to, msg_nto);
        }

        cursor->tcur = tnext;
        iproc_message_iter_select(it, 0);
        has_next = 1;
    }

    return has_next;
}

int64_t
iproc_cursor_time (iproc_cursor *cursor)
{
    if (!cursor)
        return INT64_MIN;

    return cursor->tcur;
}

iproc_history *
iproc_cursor_history (iproc_cursor *cursor)
{
    if (!cursor)
        return NULL;

    return cursor->history;
}

int64_t
iproc_cursor_nmsg (iproc_cursor *cursor)
{
    if (!cursor)
        return 0;

    return iproc_message_iter_ntie(cursor->it);
}

void
iproc_cursor_select_msg (iproc_cursor *cursor,
                         int64_t       i)
{
    assert(cursor);
    assert(i >= 0);
    assert(i < iproc_cursor_nmsg(cursor));
    iproc_message_iter_select(cursor->it, i);
}

int64_t
iproc_cursor_msg_from (iproc_cursor *cursor)
{
    assert(cursor);
    return iproc_message_iter_from(cursor->it);
}

int64_t
iproc_cursor_msg_nto (iproc_cursor *cursor)
{
    assert(cursor);
    return iproc_message_iter_nto(cursor->it);
}

int64_t *
iproc_cursor_msg_to (iproc_cursor *cursor)
{
    assert(cursor);
    return iproc_message_iter_to(cursor->it);
}

