
require(RUnit)
require(iproc)
data(enron)

frame <- coef <- has.loops <- m <- NULL
msgs <- it <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <- actors( ~ -1, data.frame(matrix(NA, nrow(a), 0)))
    receivers <- senders
    receive.intervals <- 3600 * 2^seq(-6, 14)
    frame <<- iproc.frame(senders, receivers, receive.intervals = receive.intervals)


    msgs <<- messages(enron)
    it <<- cursor(msgs)

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
    while (advance(it)) {
        for (i in from(it)) {
            lw <- t(mul(frame, coef, sender = i, it))
            if (!has.loops(m)) {
                lw[i] <- -Inf
            }

            jmax <- which.max(lw)
            lw.max <- lw[jmax]
            lw1 <- lw - lw.max
            scale <- log1p(sum(exp(lw1[-jmax])))
            lp <- lw1 - scale
            
            checkEquals(log.probs(m, i, it), lp)
        }
    }
}

test.probs <- function() {
    while (advance(it)) {
        for (i in from(it)) {
            checkEquals(probs(m, i, it), exp(log.probs(m, i, it)))
        }
    }
}
