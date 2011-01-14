
messages.default <- function(time, from, to, nto = rep(1, length(from))) {
    time <- as.integer(time)
    from <- as.integer(from)
    to <- as.integer(to)
    nto <- as.integer(nto)
    .Call("Riproc_messages_new", time, from, to, nto)
}

size.messages <- function(messages) {
    .Call("Riproc_messages_size", messages)
}

time.messages <- function(messages) {
    .Call("Riproc_messages_time", messages)
}

from.messages <- function(messages) {
    .Call("Riproc_messages_from", messages)
}

to.messages <- function(messages) {
    .Call("Riproc_messages_to", messages)
}
    
nto.messages <- function(messages) {
    .Call("Riproc_messages_nto", messages)
}

max.from.messages <- function(messages) {
    .Call("Riproc_messages_max_from", messages)
}

max.to.messages <- function(messages) {
    .Call("Riproc_messages_max_to", messages)
}

max.nto.messages <- function(messages) {
    .Call("Riproc_messages_max_nto", messages)
}
