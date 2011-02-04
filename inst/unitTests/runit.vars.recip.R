
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- receive.intervals <- v <- NULL
msgs <- it <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <<- actors(matrix(numeric(), size(a), 0))
    receivers <<- senders
    receive.intervals <<- 3600 * 2^seq(-6, 14)
    v <<- vars(senders, receivers, receive.intervals = receive.intervals)

    msgs <<- messages(enron)
    it <<- cursor(msgs)

    set.seed(0)
}

.tearDown <- function() {
    senders <<- receivers <<- v <<- NULL
    msgs <<- it <<- NULL
    gc()
}

test.dim <- function() {
    checkEquals(dim(v), length(receive.intervals))
}

test.as.matrix <- function() {
    tlast <- matrix(-Inf, size(senders), size(receivers))
    while(advance(it)) {
        tcur <- time(it)
        for (i in from(it)) {
            delta <- tcur - tlast[,i]
            x <- matrix(0.0, size(receivers), dim(v))
            for (j in seq_len(size(receivers))) {
                int <- which(delta[j] <= receive.intervals)
                if (any(int)) {
                    x[j, min(int)] <- 1.0
                }
            }

            checkEquals(as.matrix(v, sender = i, it), x)
        }

        for (t in seq_len(nties(it))) {
            tlast[from(it)[t], to(it)[[t]]] <- tcur
        }
    }
}


test.mul <- function() {
    x <- sample(-2:2, dim(v), replace = TRUE)
    while (advance(it)) {
        for (i in from(it)) {
            checkEquals(mul(v, x, sender = i, it), as.matrix(v, sender = i, it) %*% x)
        }
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nreceiver(v), replace = TRUE)
    while (advance(it)) {
        for (i in from(it)) {
            checkEquals(tmul(v, x, sender = i, it), t(as.matrix(v, sender = i, it)) %*% x)
        }
    }
}
