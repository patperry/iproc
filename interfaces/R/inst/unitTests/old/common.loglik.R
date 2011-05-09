require(RUnit)
require(iproc)

data(enron)

msgs <- senders <- receivers <- recip.intervals <- design <- beta <- mdl <- it <- NULL
max.advance <- NULL

loglik.setUp <- function(senders, receivers, recip.intervals = NULL,
                         has.loops = FALSE) {
    set.seed(0)
    msgs <<- messages(time, sender.id, receiver.id, enron$messages)
    senders <<- senders
    receivers <<- receivers
    recip.intervals <<- recip.intervals
    design <<- iproc.design(senders, receivers, recip.intervals = recip.intervals)
    beta <<- sample(-2:2, ncol(design), replace = TRUE)
    mdl <<- model(design, beta, has.loops = has.loops)
    it <<- cursor(msgs)
    max.advance <<- 100
}

test.value <- function() {
    value <- 0.0
    ll <- loglik(mdl)
    n <- 0
    
    while(advance(it)) {
        n <- n + 1
        for (tie in seq_len(nrow(it))) {
            msg.from <- from(it)[[tie]]
            msg.to <- to(it)[[tie]]

            lp <- as.vector(log.probs(mdl, msg.from, it))
            value <- value + sum(lp[msg.to])
        }

        insert(ll, it)
        checkEquals(value(ll), value)

        if (n == max.advance)
            break
    }
}

test.grad <- function() {
    ll <- loglik(mdl)
    mdl0 <- model(design, beta, has.loops = TRUE)
    grad <- rep(0.0, ncol(design))
    nrecv <- nrow(design)
    

    err <- c()
    
    # advance(it); advance(it); advance(it); advance(it);
    # advance(it); advance(it); advance(it); advance(it);
    # advance(it); advance(it); advance(it); advance(it);
    # advance(it);

    n <- 0
    while(advance(it)) {
        n <- n + 1
        for (tie in seq_len(nrow(it))) {
            msg.from <- from(it)[[tie]]
            msg.to <- to(it)[[tie]]
            msg.nto <- length(msg.to)

            n.expected <- msg.nto * as.vector(probs(mdl, msg.from, it))

            n.actual <- rep(0, nrecv)
            for (t in seq_along(msg.to)) {
                n.actual[msg.to[t]] <- n.actual[msg.to[t]] + 1.0
            }
            dgrad <- tmul(design, n.actual - n.expected, sender = msg.from, it)


            # w0 <- exp(mul(design, beta, msg.from))
            # w <- exp(mul(design, beta, msg.from, it)); w[msg.from] <- 0
            # suminvwt <- msg.nto * (sum(w0)/sum(w))
            # p0 <- as.vector(probs(mdl0, msg.from))
            # p <- as.vector(probs(mdl, msg.from, it))
            # p.active <- rep(0, nrecv)            
            # p.active[w != w0] <- p[w != w0]
            #
            # dp <- p.active
            # dp[w != w0] <- (msg.nto * p.active - suminvwt * p0)[w != w0]
            #
            # x <- tmul(design, n.actual, msg.from, it)
            # e1 <- suminvwt * tmul(design, p0, msg.from)
            # e2 <- tmul(design, dp, msg.from)
            # e3 <- msg.nto * (tmul(design, p, msg.from, it) - tmul(design, p, msg.from))
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
        if (n == max.advance)
            break
    }
}
