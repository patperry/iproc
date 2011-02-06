
require(RUnit)
require(iproc)

data(enron)

time <- from <- to <- msgs <- NULL

.setUp <- function() {
    time <<- as.integer(enron$messages$time)
    from <<- enron$messages$sender.id
    to <<- enron$messages$receiver.id
    msgs <<- messages(time, sender.id, receiver.id, enron$messages)
}

test.size <- function() {
    checkEquals(size(msgs), length(time))
}

test.max.from <- function() {
    checkEquals(max.from(msgs), max(from))
}

test.max.to <- function() {
    checkEquals(max.to(msgs), max(sapply(to, max)))
}

test.max.nto <- function() {
    checkEquals(max.nto(msgs), max(sapply(to, length)))
}

test.time <- function() {
    checkEquals(time(msgs), time)
}

test.from <- function() {
    checkEquals(from(msgs), from)
}

test.to <- function() {
    # the check goes a lot faster if we concatenate
    checkEquals(c(to(msgs), recursive = TRUE),
                c(to, recursive = TRUE))
}
    
test.nto <- function() {
    checkEquals(nto(msgs), sapply(to, length))
}

