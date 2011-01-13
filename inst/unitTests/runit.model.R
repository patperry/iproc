
require(RUnit)
require(iproc)

vars <- coef <- has.loops <- m <- NULL

.setUp <- function() {
    senders <- actors(c(2,1,2), matrix(1:6, 2, 3))
    receivers <- actors(1:4, rbind(2,4,6,8) %*% cbind(1,-1))
    vars <<- vars(senders, receivers)
    coef <<- c(0, 3, -2, 0, 1, -1)
    has.loops <<- FALSE
    m <<- model(vars, coef, has.loops)
}

.tearDown <- function() {
    vars <<- coef <<- has.loops <<- m <<- NULL
    gc()
}


test.vars <- function() {
    checkIdentical(vars(m), vars)
}

test.coef <- function() {
    checkEquals(coef(m), coef)
}

test.has.loops <- function() {
    checkEquals(has.loops(m), has.loops)
}

test.dim <- function() {
    checkEquals(dim(m), dim(vars))
}

test.nsender <- function() {
    checkEquals(nsender(m), nsender(vars))
}

test.nreceiver <- function() {
    checkEquals(nreceiver(m), nreceiver(vars))
}
