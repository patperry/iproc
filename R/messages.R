
messages <- function(time, from, to, data = NULL) {
    if (missing(time))
        stop("Must have a time argument")
    if (missing(from))
        stop("Must have a from argument")
    if (missing(to))
        stop("Must have a to argument")

    time <- eval(call("with", substitute(data), substitute(time)))
    from <- eval(call("with", substitute(data), substitute(from)))
    to <- eval(call("with", substitute(data), substitute(to)))    

    if (length(time) != length(from) || length(time) != length(to))
        stop("Time, from, and to must have same lengths.")

    rawtime <- unclass(time)
    if (!is.numeric(rawtime))
        stop("Time must inherit from numeric or POSIXct")
    from <- as.integer(from)
    to <- lapply(as.list(to), as.integer)

    msgs <- .Call("Riproc_messages_new", rawtime, from, to)
    attr(msgs, "attributes.time") <- attributes(time)
    msgs
}

dim.messages <- function(x, ...) {
    messages <- x
    n <- .Call("Riproc_messages_size", messages)
    p <- 3L
    c(n,p)
}

time.messages <- function(x, ...) {
    messages <- x
    time <- .Call("Riproc_messages_time", messages)
    attributes(time) <- attr(messages, "attributes.time")
    time
}

from.messages <- function(x, ...) {
    messages <- x
    .Call("Riproc_messages_from", messages)
}

to.messages <- function(x, ...) {
    messages <- x
    .Call("Riproc_messages_to", messages)
}
    
# nto.messages <- function(messages) {
#     .Call("Riproc_messages_nto", messages)
# }

# max.from.messages <- function(messages) {
#     .Call("Riproc_messages_max_from", messages)
# }

# max.to.messages <- function(messages) {
#     .Call("Riproc_messages_max_to", messages)
# }

# max.nto.messages <- function(messages) {
#     .Call("Riproc_messages_max_nto", messages)
# }
