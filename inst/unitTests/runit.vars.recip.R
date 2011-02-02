
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- receive.intervals <- v <- NULL
msgs <- it <- NULL

.setUp <- function() {
    a <- actors(enron)
    senders <<- actors(group(a), matrix(numeric(0), ngroup(a), 0))
    receivers <<- senders
    receive.intervals <<- 10^seq_len(6)
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
