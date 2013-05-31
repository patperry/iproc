
#define FOREACH_ACTOR(k, i, msg) \
	for ((k) = 0; (k) < 1 && ((i) = (msg)->from, 1); (k)++)

#include "v-itot-source.h"

const struct tvar_type *VAR_ISENDTOT = &VAR_ITOT_REP;

