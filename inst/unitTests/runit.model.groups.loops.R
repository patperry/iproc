
require(RUnit)
require(iproc)
data(enron)

vars <- coef <- has.loops <- m <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <- a
    receivers <- a
    vars <<- vars(senders, receivers)
    set.seed(0)
    coef <<- sample(-2:2, dim(vars), replace = TRUE)
    has.loops <<- TRUE
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

test.log.probs <- function() {
    for (i in seq_len(nsender(vars))) {
        lw <- t(mul(vars, coef, sender = i))
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
    for (i in seq_len(nsender(vars))) {
        checkEquals(probs(m, i), exp(log.probs(m, i)))
    }
}
