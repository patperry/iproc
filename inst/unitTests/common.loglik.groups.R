
source(file.path(runit.dir(), "common.loglik.R"), TRUE)

loglik.groups.setUp <- function(has.loops) {
    actrs <- actors(enron)
    loglik.setUp(actrs, actrs, has.loops = has.loops)
}
