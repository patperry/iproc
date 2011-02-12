# actors methods
group        <- function(object, ...) UseMethod("group")
ngroup       <- function(object, ...) UseMethod("ngroup")
group.traits <- function(object, ...) UseMethod("group.traits")

# messages methods
from         <- function(x, ...) UseMethod("from")
to           <- function(x, ...) UseMethod("to")

# iproc.design methods
iproc.design <- function(object, ...) UseMethod("iproc.design")
receivers    <- function(object, ...) UseMethod("receivers")
senders      <- function(object, ...) UseMethod("senders")

# matrix multiplicaiton
mul          <- function(x, y, ...) UseMethod("mul")
tmul         <- function(x, y, ...) UseMethod("tmul")

# cursor
advance      <- function(cursor, ...) UseMethod("advance")
cursor       <- function(object, ...) UseMethod("cursor")
finished     <- function(cursor, ...) UseMethod("finished")
reset        <- function(cursor, ...) UseMethod("reset")
started      <- function(cursor, ...) UseMethod("started")



grad         <- function(object, ...) UseMethod("grad")
has.loops    <- function(object, ...) UseMethod("has.loops")
insert       <- function(object, ...) UseMethod("insert")
loglik       <- function(object, ...) UseMethod("loglik")
log.probs    <- function(object, ...) UseMethod("log.probs")
model        <- function(object, ...) UseMethod("model")
nreceiver    <- function(object, ...) UseMethod("nreceiver")
nsender      <- function(object, ...) UseMethod("nsender")
probs        <- function(object, ...) UseMethod("probs")
value        <- function(object, ...) UseMethod("value")
