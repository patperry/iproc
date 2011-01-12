#ifndef _RIPROC_UTILS_H
#define _RIPROC_UTILS_H

#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

void *   Riproc_sexp2ptr  (SEXP     s,
                           Rboolean null_ok,
                           SEXP     tag,
                           char    *type);

Rboolean Riproc_sexpisptr (SEXP     s,
                           SEXP     tag);



#endif /* _RIPROC_UTILS_H */
