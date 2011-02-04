
iproc.frame.default <- function(senders, receivers, data, data.receivers = data,
                                receive.intervals = NULL, data.senders = data, ...)
{
    if (!inherits(senders, "actors")) {
        if (missing(data) && missing(data.senders)) {
            senders <- actors(senders)
        } else {
            senders <- actors(senders, data.senders)
        }
    }
    
    if (!inherits(receivers, "actors")) {
        if (missing(data) && missing(data.receivers)) {
            receivers <- actors(receivers)
        } else {
            receivers <- actors(receivers, data.receivers)
        }
    }

    receive.intervals.orig <- as.integer(receive.intervals)
    if (any(is.na(receive.intervals.orig) | (receive.intervals.orig <= 0)))
        stop("'receive.intervals' contains a nonpositive or NA value")
    
    receive.intervals <- unique(sort(receive.intervals.orig))
    if (!identical(receive.intervals, receive.intervals.orig))
        stop("'receive.intervals' not unique and sorted in ascending order")

    frame <- .Call("Riproc_frame_new", senders, receivers, receive.intervals)
    attr(frame, "class") <- "iproc.frame"
    attr(frame, "attributes.senders") <- attributes(senders)
    attr(frame, "attributes.receivers") <- attributes(receivers)
    frame
}

dim.iproc.frame <- function(object, ...) {
    frame <- object
    n <- .Call("Riproc_frame_nreceiver", frame)
    p <- .Call("Riproc_frame_dim", frame)
    c(n,p)
}

senders.iproc.frame <- function(object, ...) {
    frame <- object
    senders <- .Call("Riproc_frame_senders", frame)
    attributes(senders) <- attr(frame, "attributes.senders")
    senders
}

receivers.iproc.frame <- function(object, ...) {
    frame <- object
    receivers <- .Call("Riproc_frame_receivers", frame)
    attributes(receivers) <- attr(frame, "attributes.receivers")
    receivers
}

as.matrix.iproc.frame <- function(object, sender, cursor = NULL, ...) {
    frame <- object
    p <- ncol(frame)
    mul(frame, diag(p), sender, cursor)
}

mul.iproc.frame <- function(x, y, sender, cursor = NULL, ...) {
    frame <- x
    matrix <- as.matrix(y)
    storage.mode(matrix) <- "numeric"
    sender <- as.integer(sender)
    .Call("Riproc_frame_mul", frame, matrix, sender, cursor)
}

tmul.iproc.frame <- function(x, y, sender, cursor = NULL, ...) {
    frame <- x
    matrix <- as.matrix(y)
    storage.mode(matrix) <- "numeric"
    sender <- as.integer(sender)
    .Call("Riproc_frame_tmul", frame, matrix, sender, cursor)
}
