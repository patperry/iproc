require(RUnit)
require(iproc)

data(enron)

msgs <- actrs <- vrs <- beta <- mdl <- it <- NULL

.setUp <- function() {
    set.seed(0)
    msgs <<- messages(enron)
    actrs <<- actors(enron)
    vrs <<- vars(actrs)
    beta <<- sample(-2:2, dim(vrs), replace = TRUE)
    mdl <<- model(vrs, beta)
    it <<- cursor(msgs)
}

test.value <- function() {
    value <- 0.0
    ll <- loglik(mdl, msgs)
    
    while(advance(it)) {
        for (tie in seq_len(nties(it))) {
            msg.from <- from(it)[[tie]]
            msg.to <- to(it)[[tie]]

            lp <- as.vector(log.probs(mdl, msg.from, it))
            value <- value + sum(lp[msg.to])
        }
    }

    checkEquals(value(ll), value)
}

test.grad <- function() {
    ll <- loglik(mdl)
    mdl0 <- model(vrs, beta, has.loops=TRUE)
    grad <- rep(0.0, dim(vrs))
    nrecv <- nreceiver(vrs)

    err <- c()
    
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
            grad <- grad + tmul(vrs, n.actual - n.expected, sender = msg.from, it)
            
            if (FALSE) { #internally, this is how the computation is done:
                w0 <- exp(mul(vrs, beta, msg.from))
                w <- exp(mul(vrs, beta, msg.from, it))
                w[msg.from] <- 0
                p0 <- as.vector(probs(mdl0, msg.from))
                p <- as.vector(probs(mdl, msg.from, it))
                e0 <- tmul(vrs, p0, msg.from)
                dp.active <- rep(0, nrecv)
                dp.active[msg.from] <- (-p0[msg.from]) * (sum(w0)/sum(w))
                e <- msg.nto * ((sum(w0)/sum(w))*e0
                                + tmul(vrs, dp.active, msg.from))
                grad <- grad + (tmul(vrs, n.actual, msg.from, it) - e)
            }
        }
        grad <- as.vector(grad)

        insert(ll, it)
        err.new <- mean(abs(grad(ll) - grad) / (abs(grad) + 0.1))
        err <- c(err, err.new)
        
        checkTrue(err.new < 1e-10)
    }
}
