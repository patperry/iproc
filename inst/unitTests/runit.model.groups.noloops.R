
require(RUnit)
require(iproc)
data(enron)

frame <- coef <- has.loops <- m <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <- a
    receivers <- a
    frame <<- iproc.frame(senders, receivers)
    set.seed(0)
    coef <<- sample(-2:2, ncol(frame), replace = TRUE)
    has.loops <<- FALSE
    m <<- model(frame, coef, has.loops)
}

.tearDown <- function() {
    frame <<- coef <<- has.loops <<- m <<- NULL
    gc()
}


test.iproc.frame <- function() {
    checkIdentical(iproc.frame(m), frame)
}

test.coef <- function() {
    checkEquals(coef(m), coef)
}

test.has.loops <- function() {
    checkEquals(has.loops(m), has.loops)
}

test.dim <- function() {
    checkEquals(dim(m), ncol(frame))
}

test.nsender <- function() {
    checkEquals(nsender(m), nrow(senders(frame)))
}

test.nreceiver <- function() {
    checkEquals(nreceiver(m), nrow(frame))
}

test.log.probs <- function() {
    for (i in seq_len(nrow(senders(frame)))) {
        lw <- t(mul(frame, coef, sender = i))
        if (!has.loops(m)) {
            lw[i] <- -Inf
        }

        jmax <- which.max(lw)
        lw.max <- lw[jmax]
        lw1 <- lw - lw.max
        scale <- log1p(sum(exp(lw1[-jmax])))
        lp <- lw1 - scale

        checkEquals(log.probs(m, i), lp)
    }
}

test.probs <- function() {
    for (i in seq_len(nrow(senders(frame)))) {
        checkEquals(probs(m, i), exp(log.probs(m, i)))
    }
}
