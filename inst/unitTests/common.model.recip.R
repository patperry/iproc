
source(file.path(runit.dir(), "common.model.R"), TRUE)

model.recip.setUp <- function(has.loops) {
    actrs <- actors(~ gender *seniority * department, enron$employees)
    actrs0 <- actors( ~ -1, data.frame(matrix(NA, nrow(actrs), 0)))
    receive.intervals <- 3600 * 2^seq(-6, 14)
    model.setUp(actrs0, actrs0, receive.intervals = receive.intervals,
                has.loops = has.loops)
}
