
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- v <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <<- a
    groups <- c(1, 1, 2, 2, 1, 2)
    group.traits <- rbind(1, -1) %*% cbind(2, 4)
    receivers <<- actors(group.traits[groups,,drop=FALSE])
    set.seed(0)
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
    checkEquals(nsender(v), nrow(senders))
}

test.nreceiver <- function() {
    checkEquals(nreceiver(v), nrow(receivers))
}

test.dim <- function() {
    checkEquals(dim(v), ncol(senders) * ncol(receivers))
}

test.as.matrix <- function() {
    s <- as.matrix(senders)
    r <- as.matrix(receivers)
    for (i in seq_len(nsender(v))) {
        checkEquals(as.matrix(v, sender = i), kronecker(r, s[i,,drop = FALSE]))
    }
}

test.mul <- function() {
    x <- sample(-2:2, dim(v), replace = TRUE)
    for (i in seq_len(nsender(v))) {
        checkEquals(mul(v, x, sender = i), as.matrix(v, sender = i) %*% x)
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nreceiver(v), replace = TRUE)
    for (i in seq_len(nsender(v))) {
        checkEquals(tmul(v, x, sender = i), t(as.matrix(v, sender = i)) %*% x)
    }
}
