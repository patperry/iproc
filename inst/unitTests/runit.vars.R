
require(RUnit)
require(iproc)

senders <- receivers <- v <- NULL

.setUp <- function() {
    senders <<- actors(c(2,1,2), matrix(1:6, 2, 3))
    receivers <<- actors(1:4, rbind(2,4,6,8) %*% cbind(1,-1))
    v <<- vars(senders, receivers)
}

.tearDown <- function() {
    senders <<- receivers <<- v <<- NULL
    gc()
}


test.senders <- function() {
    checkIdentical(senders(v), senders)
}

test.receivers <- function() {
    checkIdentical(receivers(v), receivers)
}

test.nsender <- function() {
    checkEquals(nsender(v), size(senders))
}

test.nreceiver <- function() {
    checkEquals(nreceiver(v), size(receivers))
}

test.dim <- function() {
    checkEquals(dim(v), dim(senders) * dim(receivers))
}
