
#define HISTORY_ACTOR1(history, i) \
	history_send((history), (i))

#define FOREACH_ACTOR1(it, h, msg) \
	for ((it) = 0; \
	     (it) < (msg)->nto \
	     && ((h) = (msg)->to[it], 1); \
	     (it)++)

#define MSG_WEIGHT1(msg) (1.0/(msg)->nto)


#define HISTORY_ACTOR2(history, h) \
	history_send((history), (h))

#define FOREACH_ACTOR2(it, j, msg) \
	for ((it) = 0; \
	     (it) < (msg)->nto \
	     && ((j) = (msg)->to[it], 1); \
	     (it)++)

#define MSG_WEIGHT2(msg) (1.0/(msg)->nto)


#include "v2-ntriad-source.h"

const struct tvar2_type *VAR2_NSEND2 = &VAR2_NTRIAD_REP;

