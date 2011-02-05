
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- receive.intervals <- f <- NULL
msgs <- it <- NULL
max.advance <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <<- actors(~ -1, data.frame(matrix(NA, nrow(a), 0)))
    receivers <<- senders
    receive.intervals <<- 3600 * 2^seq(-6, 14)
    f <<- iproc.frame(senders, receivers, receive.intervals = receive.intervals)

    msgs <<- messages(enron)
    it <<- cursor(msgs)

    max.advance <<- 100
    
    set.seed(0)
}

.tearDown <- function() {
    senders <<- receivers <<- f <<- NULL
    msgs <<- it <<- NULL
    gc()
}

test.dim <- function() {
    checkEquals(ncol(f), length(receive.intervals))
}

test.as.matrix <- function() {
    tlast <- matrix(-Inf, nrow(senders), nrow(receivers))
    n <- 0
    while(advance(it)) {
        n <- n + 1
        tcur <- time(it)
        for (i in from(it)) {
            delta <- tcur - tlast[,i]
            x <- matrix(0.0, nrow(receivers), ncol(f))
            for (j in seq_len(nrow(receivers))) {
                int <- which(delta[j] <= receive.intervals)
                if (any(int)) {
                    x[j, min(int)] <- 1.0
                }
            }

            checkEquals(as.matrix(f, sender = i, it), x)
        }

        for (t in seq_len(nties(it))) {
            tlast[from(it)[t], to(it)[[t]]] <- tcur
        }

        if (n == max.advance)
            break
    }
}


test.mul <- function() {
    x <- sample(-2:2, ncol(f), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(mul(f, x, sender = i, it),
                        as.matrix(f, sender = i, it) %*% x)
        }
        if (n == max.advance)
            break
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nrow(f), replace = TRUE)
    n <- 0
    while (advance(it)) {
        n <- n + 1
        for (i in from(it)) {
            checkEquals(tmul(f, x, sender = i, it),
                        t(as.matrix(f, sender = i, it)) %*% x)
        }

        if (n == max.advance)
            break
    }
}
