#ifndef _RIPROC_FIT_H
#define _RIPROC_FIT_H


#include <R.h>
#include <Rdefines.h>


/* Call once to initialize library */
void Riproc_fit_init (DllInfo *info);


SEXP Riproc_fit (SEXP Rmodel0,
                 SEXP Rmessages,
                 SEXP Rpenalty);

#endif /* _RIPROC_FIT_H */
