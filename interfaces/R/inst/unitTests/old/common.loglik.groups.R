
source(file.path(runit.dir(), "common.loglik.R"), TRUE)

loglik.groups.setUp <- function(has.loops) {
    actrs <- actors(~ gender *seniority * department, enron$employees)
    loglik.setUp(actrs, actrs, has.loops = has.loops)
}
