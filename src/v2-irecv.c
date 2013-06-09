
#define HISTORY_ACTOR(h, i) \
	history_recv((h), (i))

#define FOREACH_ACTOR(it, j, msg) \
	for ((it) = 0; \
	     (it) < 1 \
	     && ((j) = (msg)->from, 1); \
	     (it)++)

#include "v2-idyad-source.h"

const struct tvar2_type *VAR2_IRECV = &VAR2_IDYAD_REP;

