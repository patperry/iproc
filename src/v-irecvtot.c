
#define FOREACH_ACTOR(k, i, msg) \
	for ((k) = 0; (k) < (msg)->nto && ((i) = (msg)->to[k], 1); (k)++)

#include "v-itot-source.h"

const struct tvar_type *VAR_IRECVTOT = &VAR_ITOT_REP;
