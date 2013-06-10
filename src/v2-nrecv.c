
#define HISTORY_ACTOR(h, i) \
	history_recv((h), (i))

#define FOREACH_ACTOR(it, j, msg) \
	for ((it) = 0; \
	     (it) < 1 \
	     && ((j) = (msg)->from, 1); \
	     (it)++)

#define MSG_WEIGHT(msg) (1.0/(msg)->nto)

#include "v2-ndyad-source.h"

const struct tvar2_type *VAR2_NRECV = &VAR2_NDYAD_REP;

