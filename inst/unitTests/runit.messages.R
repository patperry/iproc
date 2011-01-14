
require(RUnit)
require(iproc)

data(enron)

time <- from <- to <- msgs <- NULL

.setUp <- function() {
    n <- nrow(enron$messages)
    time <<- enron$messages$unix.time[1:n]
    from <<- enron$messages$sender.id[1:n]
    to <<- enron$messages$receiver.id[1:n]
    msgs <<- messages(time, from, to)
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
    checkEquals(to(msgs), to)
}
    
test.nto <- function() {
    checkEquals(nto(msgs), sapply(to, length))
}

