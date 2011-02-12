
require(RUnit)
require(iproc)
data(enron)

design <- coef <- has.loops <- m <- NULL
msgs <- it <- NULL
max.advance <- NULL

model.setUp <- function(senders, receivers,
                        receive.intervals = NULL,
                        has.loops = FALSE) {
    design <<- iproc.design(senders, receivers, receive.intervals = receive.intervals)
    set.seed(0)
    coef <<- sample(-2:2, ncol(design), replace = TRUE)
    has.loops <<- has.loops
    m <<- model(design, coef, has.loops)
    msgs <<- messages(time, sender.id, receiver.id, enron$messages)
    it <<- cursor(msgs)
    max.advance <<- 100
}

.tearDown <- function() {
    design <<- coef <<- has.loops <<- m <<- NULL
    gc()
}

test.iproc.design <- function() {
    checkIdentical(iproc.design(m), design)
}

test.coef <- function() {
    checkEquals(coef(m), coef)
}

test.has.loops <- function() {
    checkEquals(has.loops(m), has.loops)
}

test.dim <- function() {
    checkEquals(dim(m), ncol(design))
}

test.nsender <- function() {
    checkEquals(nsender(m), nrow(senders(design)))
}

test.nreceiver <- function() {
    checkEquals(nreceiver(m), nrow(design))
}

test.log.probs0 <- function() {
    for (i in seq_len(nrow(senders(design)))) {
        lw <- t(mul(design, coef, sender = i))
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

test.log.probs <- function() {
    n <- 0
    while (advance(it)) {
        n <- n + 1

        for (tie in seq_len(nrow(it))) {
            i <- from(it)[[tie]]

            lw <- t(mul(design, coef, sender = i, it))
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

        if (n == max.advance)
            break
    }
}

test.probs0 <- function() {
    for (i in seq_len(nrow(senders(design)))) {
        checkEquals(probs(m, i), exp(log.probs(m, i)))
    }
}

test.probs <- function() {
    n <- 0
    while (advance(it)) {
        n <- n + 1

        for (tie in seq_len(nrow(it))) {
            i <- from(it)[[tie]]
            checkEquals(probs(m, i), exp(log.probs(m, i)))
        }
        
        if (n == max.advance)
            break
    }
}
