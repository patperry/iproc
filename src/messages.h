#ifndef _IPROC_MESSAGES_H
#define _IPROC_MESSAGES_H

#include <stdbool.h>
#include <stdint.h>
#include "darray.h"
#include "history.h"
#include "refcount.h"

typedef struct _iproc_message      iproc_message;
typedef struct _iproc_messages     iproc_messages;
typedef struct _iproc_message_iter iproc_message_iter;

struct _iproc_message {
    double   time;
    int64_t  from;
    int64_t  ito;
    int64_t  nto;
};

struct _iproc_messages {
    double tcur;
    struct darray  array;
    struct darray  recipients;
    iproc_refcount refcount;
    int64_t max_from;
    int64_t max_to;
    int64_t max_nto;
};

struct _iproc_message_iter {
    iproc_messages *messages;
    iproc_history  *history;
    int64_t         offset;
    int64_t         ntie;
    iproc_message  *message;
    bool            finished;
    iproc_refcount  refcount;
};

iproc_messages * iproc_messages_new        ();
iproc_messages * iproc_messages_ref        (iproc_messages *msgs);
void             iproc_messages_unref      (iproc_messages *msgs);

int64_t          iproc_messages_size       (iproc_messages *msgs);
void             iproc_messages_advance_to (iproc_messages *msgs,
                                            double          t);


void             iproc_messages_insert     (iproc_messages *msgs,
                                            int64_t         from,
                                            int64_t         to);
void             iproc_messages_insertm    (iproc_messages *msgs,
                                            int64_t         from,
                                            int64_t        *to,
                                            int64_t         nto);

int64_t          iproc_messages_max_from   (iproc_messages *msgs);
int64_t          iproc_messages_max_to     (iproc_messages *msgs);
int64_t          iproc_messages_max_nto    (iproc_messages *msgs);


iproc_message_iter * iproc_message_iter_new    (iproc_messages     *msgs);
iproc_message_iter * iproc_message_iter_ref    (iproc_message_iter *it);
void                 iproc_message_iter_unref  (iproc_message_iter *it);

double               iproc_message_iter_time   (iproc_message_iter *it);
int64_t              iproc_message_iter_ntie   (iproc_message_iter *it);
void                 iproc_message_iter_select (iproc_message_iter *it,
                                                int64_t             tie);

int64_t              iproc_message_iter_from   (iproc_message_iter *it);
int64_t              iproc_message_iter_nto    (iproc_message_iter *it);
int64_t *            iproc_message_iter_to     (iproc_message_iter *it);

iproc_history *      iproc_message_iter_history (iproc_message_iter *it);

void                 iproc_message_iter_reset    (iproc_message_iter *it);
bool                 iproc_message_iter_next     (iproc_message_iter *it);
bool                 iproc_message_iter_started  (iproc_message_iter *it);
bool                 iproc_message_iter_finished (iproc_message_iter *it);

#endif /* _IPROC_MESSAGES_H */
