require(RUnit)
require(iproc)

data(enron)

msgs <- actrs <- vrs <- beta <- has.loops <- mdl <- it <- NULL

.setUp <- function() {
    set.seed(0)
    msgs <<- messages(enron)
    a <- actors(enron)
    actrs <<- actors(group(a), matrix(numeric(0), ngroup(a), 0))
    receive.intervals <- 3600 * 2^seq(-6, 14)
    vrs <<- vars(actrs, actrs, receive.intervals = receive.intervals)
    beta <<- sample(-2:2, dim(vrs), replace = TRUE)
    has.loops <- FALSE
    
    mdl <<- model(vrs, beta, has.loops)
    it <<- cursor(msgs)
}

test.value <- function() {
    value <- 0.0
    ll <- loglik(mdl)

    checkEquals(value(ll), value)
    
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
}

test.grad <- function() {
    ll <- loglik(mdl)
    grad <- rep(0.0, dim(vrs))
    nrecv <- nreceiver(vrs)

    checkEquals(grad(ll), grad)
    
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
            dgrad <- tmul(vrs, n.actual - n.expected, sender = msg.from, it)
            grad <- grad + dgrad

        }
        grad <- as.vector(grad)

        insert(ll, it)
        err.new <- mean(abs(grad(ll) - grad) / (abs(grad) + 0.1))

        checkTrue(err.new < 1e-10)
    }
}
