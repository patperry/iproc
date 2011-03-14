
source(file.path(runit.dir(), "common.model.R"), TRUE)

model.all.setUp <- function(has.loops) {
    actrs <- actors(~ gender *seniority * department, enron$employees)
    recip.intervals <- 3600 * 2^seq(-6, 14)
    model.setUp(actrs, actrs, recip.intervals = recip.intervals,
                has.loops = has.loops)
}
