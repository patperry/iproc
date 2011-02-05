
source(file.path(runit.dir(), "common.model.R"), TRUE)

model.groups.setUp <- function(has.loops) {
    a <- actors(enron)
    senders <- a
    receivers <- a
    model.setUp(senders, receivers, has.loops = has.loops)
}
