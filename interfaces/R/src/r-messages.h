#ifndef _RIPROC_MESSAGES_H
#define _RIPROC_MESSAGES_H

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "messages.h"

/* Call once to initialize library */
void Riproc_messages_init(DllInfo * info);

/* External functions */
SEXP Riproc_messages_new(SEXP Rtime, SEXP Rfrom, SEXP Rto);
SEXP Riproc_messages_size(SEXP Rmsgs);
SEXP Riproc_messages_time(SEXP Rmsgs);
SEXP Riproc_messages_from(SEXP Rmsgs);
SEXP Riproc_messages_to(SEXP Rmsgs);
SEXP Riproc_messages_nto(SEXP Rmsgs);

SEXP Riproc_messages_max_from(SEXP Rmsgs);
SEXP Riproc_messages_max_to(SEXP Rmsgs);
SEXP Riproc_messages_max_nto(SEXP Rmsgs);

/* Internal use only */
struct messages *Riproc_to_messages(SEXP Rmsgs);
SEXP Riproc_from_messages(struct messages * msgs);

#endif /* _RIPROC_MESSAGES_H */
