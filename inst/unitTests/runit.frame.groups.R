
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- f <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <<- a
    groups <- c(1, 1, 2, 2, 1, 2)
    group.traits <- rbind(1, -1) %*% cbind(2, 4)
    receivers <<- actors(~ . -1, data.frame(group.traits[groups,,drop=FALSE]))
    set.seed(0)
    f <<- iproc.frame(senders, receivers)
}

.tearDown <- function() {
    senders <<- receivers <<- f <<- NULL
    gc()
}


test.senders <- function() {
    checkIdentical(senders(f), senders)
}

test.receivers <- function() {
    checkIdentical(receivers(f), receivers)
}

test.dim <- function() {
    checkEquals(nrow(f), nrow(receivers))
    checkEquals(ncol(f), ncol(senders) * ncol(receivers))
}

test.as.matrix <- function() {
    s <- as.matrix(senders)
    r <- as.matrix(receivers)
    for (i in seq_len(nrow(senders(f)))) {
        checkEquals(as.matrix(f, sender = i),
                    kronecker(r, s[i,,drop = FALSE]))
    }
}

test.mul <- function() {
    x <- sample(-2:2, ncol(f), replace = TRUE)
    for (i in seq_len(nrow(senders(f)))) {
        checkEquals(mul(f, x, sender = i),
                    as.matrix(f, sender = i) %*% x)
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nrow(f), replace = TRUE)
    for (i in seq_len(nrow(senders(f)))) {
        checkEquals(tmul(f, x, sender = i),
                    t(as.matrix(f, sender = i)) %*% x)
    }
}
