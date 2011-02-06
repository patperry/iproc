
require(RUnit)
require(iproc)

data(enron)

msgs <<- it <<- NULL


.setUp <- function() {
    msgs <<- messages(time, sender.id, receiver.id, enron$messages[1:100,])
    it <<- cursor(msgs)
}

test.started <- function() {
    checkEquals(started(it), FALSE)
    
    advance(it)
    checkEquals(started(it), TRUE)
    
    while (advance(it)) { }
    checkEquals(started(it), TRUE)

    advance(it)
    checkEquals(started(it), TRUE)
}

test.finished <- function() {
    checkEquals(finished(it), FALSE)

    if (advance(it))
        checkEquals(finished(it), FALSE)

    while (advance(it)) { }
    checkEquals(finished(it), TRUE)

    advance(it)
    checkEquals(finished(it), TRUE)
}

test.advance <- function() {
    n <- 0
    while (advance(it)) {
        n <- n + nties(it)
    }
    checkEquals(n, size(msgs))
}
