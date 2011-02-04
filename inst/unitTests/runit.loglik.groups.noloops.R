require(RUnit)
require(iproc)

data(enron)

msgs <- actrs <- frame <- beta <- has.loops <- mdl <- it <- NULL

.setUp <- function() {
    set.seed(0)
    msgs <<- messages(enron)
    actrs <<- actors(enron)
    frame <<- iproc.frame(actrs, actrs)
    beta <<- sample(-2:2, ncol(frame), replace = TRUE)
    has.loops <- FALSE
    
    mdl <<- model(frame, beta, has.loops)
    it <<- cursor(msgs)
}

test.value <- function() {
    value <- 0.0
    ll <- loglik(mdl)
    
    while(advance(it)) {
        for (tie in seq_len(nties(it))) {
            msg.from <- from(it)[[tie]]
            msg.to <- to(it)[[tie]]

            lp <- as.vector(log.probs(mdl, msg.from, it))
            value <- value + sum(lp[msg.to])
        }

        insert(ll, it)
        checkEquals(value(ll), value)
    }

    checkEquals(value(ll), value)
}

test.grad <- function() {
    ll <- loglik(mdl)
    grad <- rep(0.0, ncol(frame))
    nrecv <- nrow(frame)
    
    while(advance(it)) {
        for (tie in seq_len(nties(it))) {
            msg.from <- from(it)[[tie]]
            msg.to <- to(it)[[tie]]
            msg.nto <- length(msg.to)

            n.expected <- msg.nto * as.vector(probs(mdl, msg.from, it))

            n.actual <- rep(0, nrecv)
            for (t in seq_along(msg.to)) {
                n.actual[msg.to[t]] <- n.actual[msg.to[t]] + 1.0
            }
            dgrad <- tmul(frame, n.actual - n.expected, sender = msg.from, it)
            grad <- grad + dgrad

        }
        grad <- as.vector(grad)

        insert(ll, it)
        err.new <- mean(abs(grad(ll) - grad) / (abs(grad) + 0.1))

        checkTrue(err.new < 1e-10)
    }
}
