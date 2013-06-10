
#define HISTORY_ACTOR(h, i) \
	history_send((h), (i))

#define FOREACH_ACTOR(it, j, msg) \
	for ((it) = 0; \
	     (it) < (msg)->nto \
	     && ((j) = (msg)->to[it], 1); \
	     (it)++)

#define MSG_WEIGHT(msg) (1.0/(msg)->nto)

#include "v2-ndyad-source.h"

const struct tvar2_type *VAR2_NSEND = &VAR2_NDYAD_REP;

