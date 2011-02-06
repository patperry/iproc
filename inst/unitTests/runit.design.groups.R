
require(RUnit)
require(iproc)
data(enron)

senders <- receivers <- design <- NULL

.setUp <- function() {
    a <- actors(~ gender *seniority * department, enron$employees)
    senders <<- a
    groups <- c(1, 1, 2, 2, 1, 2)
    group.traits <- rbind(1, -1) %*% cbind(2, 4)
    receivers <<- actors(~ . -1, data.frame(group.traits[groups,,drop=FALSE]))
    set.seed(0)
    design <<- iproc.design(senders, receivers)
}

.tearDown <- function() {
    senders <<- receivers <<- design <<- NULL
    gc()
}


test.senders <- function() {
    checkIdentical(senders(design), senders)
}

test.receivers <- function() {
    checkIdentical(receivers(design), receivers)
}

test.dim <- function() {
    checkEquals(nrow(design), nrow(receivers))
    checkEquals(ncol(design), ncol(senders) * ncol(receivers))
}

test.as.matrix <- function() {
    s <- as.matrix(senders)
    r <- as.matrix(receivers)
    for (i in seq_len(nrow(senders(design)))) {
        checkEquals(as.matrix(design, sender = i),
                    kronecker(r, s[i,,drop = FALSE]))
    }
}

test.mul <- function() {
    x <- sample(-2:2, ncol(design), replace = TRUE)
    for (i in seq_len(nrow(senders(design)))) {
        checkEquals(mul(design, x, sender = i),
                    as.matrix(design, sender = i) %*% x)
    }
}

test.tmul <- function() {
    x <- sample(-2:2, nrow(design), replace = TRUE)
    for (i in seq_len(nrow(senders(design)))) {
        checkEquals(tmul(design, x, sender = i),
                    t(as.matrix(design, sender = i)) %*% x)
    }
}
