
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- receive.intervals <- v <- v.group <- v.recip <- NULL
msgs <- it <- NULL

.setUp <- function() {
    a <- actors(enron)
    a0 <- actors(group(a), matrix(numeric(0), ngroup(a), 0))
    senders <<- a
    receivers <<- a
    receive.intervals <<- 3600 * 2^seq(-6, 14)
    v <<- vars(senders, receivers, receive.intervals = receive.intervals)
    v.group <<- vars(senders, receivers)
    v.recip <<- vars(a0, a0, receive.intervals = receive.intervals)

    msgs <<- messages(enron)
    it <<- cursor(msgs)

    set.seed(0)
}

test.dim <- function() {
    checkEquals(dim(v), dim(v.group) + dim(v.recip))
}

test.as.matrix <- function() {
    n <- 0
    while(advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(as.matrix(v, sender = i, it),
                        cbind(as.matrix(v.recip, sender = i, it),
                              as.matrix(v.group, sender = i, it)))
        }

        if (n == 500)
            break
    }
}


test.mul <- function() {
    x <- sample(-2:2, dim(v), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(mul(v, x, sender = i, it), as.matrix(v, sender = i, it) %*% x)
        }

        if (n == 500) break
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nreceiver(v), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <-  n + 1
        
        for (i in from(it)) {
            checkEquals(tmul(v, x, sender = i, it), t(as.matrix(v, sender = i, it)) %*% x)
        }

        if (n == 500) break
    }
}
