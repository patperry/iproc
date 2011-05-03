#ifndef _VNRECV_H
#define _VNRECV_H

#include "array.h"
#include "design.h"


/* If intvls = { delta[0], delta[1], ..., delta[K-1] }, then the K + 1 variables
 * are x{k}[t,i,j] = #{ j -> i in [t - delta[k], t - delta[k-1]) }
 * for k in [0 .. K], where delta[-1] = 0 and delta[K] = Infty.
 */

struct vnrecv {
	struct dyad_var dyad_var;
	struct array intvls;
};

bool vnrecv_init(struct vnrecv *v, const double *intvls, ssize_t n);
void vnrecv_deinit(struct vnrecv *v);

#endif /* _VNRECV_H */