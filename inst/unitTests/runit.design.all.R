
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- receive.intervals <- NULL
design <- design.group <- design.recip <- NULL
msgs <- it <- NULL
max.advance <- NULL

.setUp <- function() {
    a <- actors(~ gender *seniority * department, enron$employees)
    a0 <- actors(~ -1, data.frame(matrix(NA, nrow(a), 0)))
    senders <<- a
    receivers <<- a
    receive.intervals <<- 3600 * 2^seq(-6, 14)
    design <<- iproc.design(senders, receivers, receive.intervals = receive.intervals)
    design.group <<- iproc.design(senders, receivers)
    design.recip <<- iproc.design(a0, a0, receive.intervals = receive.intervals)

    msgs <<- messages(time, sender.id, receiver.id, enron$messages)
    it <<- cursor(msgs)

    max.advance <<- 100

    set.seed(0)
}

test.dim <- function() {
    checkEquals(ncol(design), ncol(design.group) + ncol(design.recip))
}

test.as.matrix <- function() {
    n <- 0
    while(advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(as.matrix(design, sender = i, it),
                        cbind(as.matrix(design.recip, sender = i, it),
                              as.matrix(design.group, sender = i, it)))
        }

        if (n == max.advance)
            break
    }
}


test.mul <- function() {
    x <- sample(-2:2, ncol(design), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(mul(design, x, sender = i, it),
                        as.matrix(design, sender = i, it) %*% x)
        }

        if (n == max.advance)
            break
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nrow(design), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <-  n + 1
        
        for (i in from(it)) {
            checkEquals(tmul(design, x, sender = i, it),
                        t(as.matrix(design, sender = i, it)) %*% x)
        }

        if (n == max.advance)
            break
    }
}
