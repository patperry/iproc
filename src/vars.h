#ifndef VARS_H
#define VARS_H

#include "design.h"
#include "design2.h"
#include "frame.h"
#include <stdarg.h>



/* Indicator 1{ j -> i in (-Infty, t) } */
extern const struct tvar2_type *VAR2_IRECV;

/* Indicator 1{ i -> j in (-Infty, t) } */
extern const struct tvar2_type *VAR2_ISEND;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct tvar2_type *VAR2_NRECV;

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ i -> j in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct tvar2_type *VAR2_NSEND;

extern const struct tvar2_type *VAR2_IRECV2;
extern const struct tvar2_type *VAR2_NRECV2;

extern const struct tvar2_type *VAR2_ISEND2;
extern const struct tvar2_type *VAR2_NSEND2;

extern const struct tvar2_type *VAR2_ISIB;
extern const struct tvar2_type *VAR2_NSIB;

extern const struct tvar2_type *VAR2_ICOSIB;
extern const struct tvar2_type *VAR2_NCOSIB;

char **var_name_alloc(const char *name, size_t len); /* DEPRECATED */
char **var_names_alloc(const char *name, size_t len, size_t n); /* DEPRECATED */
char **var_names_alloc2(const char *name, size_t len, size_t m, size_t n); /* DEPRECATED */
void var_names_free(char **names); /* DEPRECATED */

#endif /* VARS_H */
