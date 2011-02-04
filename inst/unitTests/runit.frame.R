
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- receive.intervals <- f <- f.group <- f.recip <- NULL
msgs <- it <- NULL

.setUp <- function() {
    a <- actors(enron)
    a0 <- actors(~ -1, data.frame(matrix(NA, nrow(a), 0)))
    senders <<- a
    receivers <<- a
    receive.intervals <<- 3600 * 2^seq(-6, 14)
    f <<- iproc.frame(senders, receivers, receive.intervals = receive.intervals)
    f.group <<- iproc.frame(senders, receivers)
    f.recip <<- iproc.frame(a0, a0, receive.intervals = receive.intervals)

    msgs <<- messages(enron)
    it <<- cursor(msgs)

    set.seed(0)
}

test.dim <- function() {
    checkEquals(ncol(f), ncol(f.group) + ncol(f.recip))
}

test.as.matrix <- function() {
    n <- 0
    while(advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(as.matrix(f, sender = i, it),
                        cbind(as.matrix(f.recip, sender = i, it),
                              as.matrix(f.group, sender = i, it)))
        }

        if (n == 500)
            break
    }
}


test.mul <- function() {
    x <- sample(-2:2, ncol(f), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(mul(f, x, sender = i, it),
                        as.matrix(f, sender = i, it) %*% x)
        }

        if (n == 500) break
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nrow(f), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <-  n + 1
        
        for (i in from(it)) {
            checkEquals(tmul(f, x, sender = i, it),
                        t(as.matrix(f, sender = i, it)) %*% x)
        }

        if (n == 500) break
    }
}
