
require(RUnit)
require(iproc)

data(enron)

time <- from <- to <- msgs <- NULL

.setUp <- function() {
    time <<- enron$messages$time
    from <<- enron$messages$sender.id
    to <<- enron$messages$receiver.id
    msgs <<- messages(time, sender.id, receiver.id, enron$messages)
}

test.nrow <- function() {
    checkIdentical(nrow(msgs), length(time))
}

# test.max.from <- function() {
#     checkIdentical(max.from(msgs), max(from))
# }

# test.max.to <- function() {
#    checkIdentical(max.to(msgs), max(sapply(to, max)))
# }

# test.max.nto <- function() {
#     checkIdentical(max.nto(msgs), max(sapply(to, length)))
# }

test.time <- function() {
    checkIdentical(time(msgs), time)
}

test.from <- function() {
    checkIdentical(from(msgs), from)
}

test.to <- function() {
    # the check goes a lot faster if we concatenate
    checkIdentical(c(to(msgs), recursive = TRUE),
                   c(to, recursive = TRUE))
}
    
test.nto <- function() {
    checkIdentical(sapply(to(msgs), length), sapply(to, length))
}
