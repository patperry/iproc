
source(file.path(runit.dir(), "common.model.R"), TRUE)

model.groups.setUp <- function(has.loops) {
    actrs <- actors(~ gender *seniority * department, enron$employees)
    senders <- actrs
    receivers <- actrs
    model.setUp(senders, receivers, has.loops = has.loops)
}
