#ifndef _VNRECV_H
#define _VNRECV_H

#include "design.h"

/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */
extern const struct var_type *VAR_TYPE_NRECV;

#endif /* _VNRECV_H */
