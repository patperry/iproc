
#define FOREACH_ACTOR(it, j, msg) \
	for ((it) = 0; \
	     (msg)->from == i \
	     && (it) < (msg)->nto \
	     && ((j) = (msg)->to[it], 1); \
	     (it)++)

#include "v2-idyad-source.h"

const struct tvar2_type *VAR2_ISEND = &VAR2_IDYAD_REP;

