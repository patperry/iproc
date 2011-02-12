
source(file.path(runit.dir(), "common.loglik.R"), TRUE)

loglik.recip.setUp <- function(has.loops) {
    actrs <- actors(~ gender *seniority * department, enron$employees)
    actrs0 <- actors(~ -1, data.frame(matrix(NA, nrow(actrs), 0)))
    recip.intervals <- 3600 * 2^seq(-6, 14)

    loglik.setUp(actrs0, actrs0, recip.intervals, has.loops = has.loops)
}
