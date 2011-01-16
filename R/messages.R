
messages.default <- function(time, from, to) {
    time <- as.integer(time)
    from <- as.integer(from)
    to <- lapply(as.list(to), as.integer)

    .Call("Riproc_messages_new", time, from, to)
}

messages.enron <- function(enron) {
    time <- as.integer(enron$messages$time)
    from <- enron$messages$sender.id
    to <- enron$messages$receiver.id
    messages.default(time, from, to)
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
