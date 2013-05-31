
#define FOREACH_ACTOR(k, i, msg) \
	for ((k) = 0; (k) < 1 && ((i) = (msg)->from, 1); (k)++)

#define MSG_WEIGHT(msg) 1.0

#include "v-ntot-source.h"

const struct tvar_type *VAR_NSENDTOT = &VAR_NTOT_REP;

