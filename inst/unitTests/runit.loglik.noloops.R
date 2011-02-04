require(RUnit)
require(iproc)

data(enron)

msgs <- actrs <- frame <- beta <- mdl <- it <- NULL

.setUp <- function() {
    set.seed(0)
    msgs <<- messages(enron)
    actrs <<- actors(enron)
    frame <<- iproc.frame(actrs, actrs,  receive.intervals = 3600 * 2^seq(-6, 14))
    beta <<- sample(-2:2, ncol(frame), replace = TRUE)
    mdl <<- model(frame, beta)
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
    mdl0 <- model(frame, beta, has.loops=TRUE)
    grad <- rep(0.0, ncol(frame))
    nrecv <- nrow(frame)
    

    err <- c()
    
    # advance(it); advance(it); advance(it); advance(it);
    # advance(it); advance(it); advance(it); advance(it);
    # advance(it); advance(it); advance(it); advance(it);
    # advance(it);
    
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


            # w0 <- exp(mul(frame, beta, msg.from))
            # w <- exp(mul(frame, beta, msg.from, it)); w[msg.from] <- 0
            # suminvwt <- msg.nto * (sum(w0)/sum(w))
            # p0 <- as.vector(probs(mdl0, msg.from))
            # p <- as.vector(probs(mdl, msg.from, it))
            # p.active <- rep(0, nrecv)            
            # p.active[w != w0] <- p[w != w0]
            #
            # dp <- p.active
            # dp[w != w0] <- (msg.nto * p.active - suminvwt * p0)[w != w0]
            #
            # x <- tmul(frame, n.actual, msg.from, it)
            # e1 <- suminvwt * tmul(frame, p0, msg.from)
            # e2 <- tmul(frame, dp, msg.from)
            # e3 <- msg.nto * (tmul(frame, p, msg.from, it) - tmul(frame, p, msg.from))
            # dgrad1 <- (((x - e1) - e2) - e3)
            #
            # grad1 <- grad + dgrad1
            
            grad <- grad + dgrad

        }
        grad <- as.vector(grad)

        insert(ll, it)
        err.new <- mean(abs(grad(ll) - grad) / (abs(grad) + 0.1))
        err <- c(err, err.new)

        #suminvwt.old <- suminvwt
        #dp.old <- dp
        #x.old <- x
        #e1.old <- e1
        #e2.old <- e2
        #e3.old <- e3
        
        checkTrue(err.new < 1e-10)

        #advance(it)
    }
}
