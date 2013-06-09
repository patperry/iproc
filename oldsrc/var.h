#ifndef VAR_H
#define VAR_H

#include <stdarg.h>
#include <stddef.h>
#include "deltaset.h"
#include "uintset.h"
#include "history.h"



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


#endif /* VAR_H */
