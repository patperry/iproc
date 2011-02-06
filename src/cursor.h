#ifndef _IPROC_CURSOR_H
#define _IPROC_CURSOR_H

#include <stdint.h>
#include "messages.h"
#include "history.h"
#include "refcount.h"


typedef struct _iproc_cursor iproc_cursor;

struct _iproc_cursor {
    double              tcur;
    iproc_history      *history;
    iproc_message_iter *it;
    iproc_refcount      refcount;
};


iproc_cursor *  iproc_cursor_new        (iproc_messages *msgs);
iproc_cursor *  iproc_cursor_ref        (iproc_cursor   *cursor);
void            iproc_cursor_unref      (iproc_cursor   *cursor);
    
void            iproc_cursor_reset      (iproc_cursor   *cursor);
int             iproc_cursor_next       (iproc_cursor   *cursor);
int             iproc_cursor_started    (iproc_cursor   *cursor);
int             iproc_cursor_finished   (iproc_cursor   *cursor);

double          iproc_cursor_time       (iproc_cursor   *cursor);
iproc_history * iproc_cursor_history    (iproc_cursor   *cursor);
int64_t         iproc_cursor_nmsg       (iproc_cursor   *cursor);

void            iproc_cursor_select_msg (iproc_cursor   *cursor,
                                         int64_t         i);
int64_t         iproc_cursor_msg_from   (iproc_cursor   *cursor);
int64_t         iproc_cursor_msg_nto    (iproc_cursor   *cursor);
int64_t *       iproc_cursor_msg_to     (iproc_cursor   *cursor);


#endif /* _IPROC_CURSOR_H */
