#ifndef _VARS_H
#define _VARS_H

#include "design.h"
#include "frame.h"
#include <stdarg.h>



/* Indicator 1{ j -> i in (-Infty, t) } */
extern const struct tvar_type *DYAD_VAR_IRECV;

/* Indicator 1{ i -> j in (-Infty, t) } */
extern const struct tvar_type *DYAD_VAR_ISEND;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct tvar_type *DYAD_VAR_NRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ i -> j in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct tvar_type *DYAD_VAR_NSEND;

//extern const struct tvar_type *RECV_VAR_IRECV2;
//extern const struct tvar_type *RECV_VAR_NRECV2;

//extern const struct tvar_type *RECV_VAR_ISEND2;
//extern const struct tvar_type *RECV_VAR_NSEND2;

//extern const struct tvar_type *RECV_VAR_ISIB;
//extern const struct tvar_type *RECV_VAR_NSIB;

//extern const struct tvar_type *RECV_VAR_ICOSIB;
//extern const struct tvar_type *RECV_VAR_NCOSIB;

char **var_name_alloc(const char *name, size_t len); /* DEPRECATED */
char **var_names_alloc(const char *name, size_t len, size_t n); /* DEPRECATED */
char **var_names_alloc2(const char *name, size_t len, size_t m, size_t n); /* DEPRECATED */
void var_names_free(char **names); /* DEPRECATED */

#endif /* _VARS_H */
