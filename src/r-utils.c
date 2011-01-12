
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "r-utils.h"

/* Adapted from Luke Tierney's Regular expressions package
 *   http://www.stat.uiowa.edu/~luke/R/regexp.html
 */
Rboolean
Riproc_sexpisptr (SEXP s,
                  SEXP tag)
{
    if (TYPEOF(s) != EXTPTRSXP || R_ExternalPtrTag(s) != tag)
        return FALSE;

    return TRUE;
}


void *
Riproc_sexp2ptr (SEXP     s,
                 Rboolean null_ok,
                 SEXP     tag,
                 char    *type)
{
    void *p;
    if (!Riproc_sexpisptr(s, tag))
        error("bad %s pointer", type);
    p = R_ExternalPtrAddr(s);
    if (!null_ok && p == NULL)
        error("null %s pointer", type);
    return p;
}
