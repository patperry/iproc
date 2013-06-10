
#define HISTORY_ACTOR1(history, i) \
	history_recv((history), (i))

#define FOREACH_ACTOR1(it, h, msg) \
	for ((it) = 0; \
	     (it) < 1 \
	     && ((h) = (msg)->from, 1); \
	     (it)++)

#define MSG_WEIGHT1(msg) (1.0/(msg)->nto)


#define HISTORY_ACTOR2(history, h) \
	history_recv((history), (h))

#define FOREACH_ACTOR2(it, j, msg) \
	for ((it) = 0; \
	     (it) < 1 \
	     && ((j) = (msg)->from, 1); \
	     (it)++)

#define MSG_WEIGHT2(msg) (1.0/(msg)->nto)


#include "v2-ntriad-source.h"

const struct tvar2_type *VAR2_NRECV2 = &VAR2_NTRIAD_REP;

