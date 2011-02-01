
require(RUnit)
require(iproc)

data(enron)

msgs <<- it <<- NULL


.setUp <- function() {
    msgs <<- messages(enron)
    it <<- cursor(msgs)
}

test.started <- function() {
    checkEquals(started(it), FALSE)
    while (advance(it)) {
        checkEquals(started(it), TRUE)
    }
    checkEquals(started(it), TRUE)
}

test.finished <- function() {
    checkEquals(finished(it), FALSE)
    while (advance(it)) {
        checkEquals(finished(it), FALSE)
    }
    checkEquals(finished(it), TRUE)
}
