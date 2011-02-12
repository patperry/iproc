
source(file.path(runit.dir(), "common.loglik.R"), TRUE)

loglik.all.setUp <- function(has.loops) {
    actrs <- actors(~ gender *seniority * department, enron$employees)
    recip.intervals <- 3600 * 2^seq(-6, 14)

    loglik.setUp(actrs, actrs, recip.intervals, has.loops = has.loops)
}
