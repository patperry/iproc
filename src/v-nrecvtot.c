
#define FOREACH_ACTOR(k, i, msg) \
	for ((k) = 0; (k) < (msg)->nto && ((i) = (msg)->to[k], 1); (k)++)

#define MSG_WEIGHT(msg) (1.0/(msg)->nto)

#include "v-ntot-source.h"

const struct tvar_type *VAR_NRECVTOT = &VAR_NTOT_REP;

