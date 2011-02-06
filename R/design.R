
iproc.design.default <- function(senders, receivers, data, data.receivers = data,
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

    receive.intervals.orig <- as.numeric(receive.intervals)
    if (any(is.na(receive.intervals.orig) | (receive.intervals.orig <= 0)))
        stop("'receive.intervals' contains a nonpositive or NA value")
    
    receive.intervals <- unique(sort(receive.intervals.orig))
    if (!identical(receive.intervals, receive.intervals.orig))
        stop("'receive.intervals' not unique and sorted in ascending order")

    design <- .Call("Riproc_design_new", senders, receivers, receive.intervals)
    attr(design, "class") <- "iproc.design"
    attr(design, "attributes.senders") <- attributes(senders)
    attr(design, "attributes.receivers") <- attributes(receivers)
    design
}

dim.iproc.design <- function(object, ...) {
    design <- object
    n <- .Call("Riproc_design_nreceiver", design)
    p <- .Call("Riproc_design_dim", design)
    c(n,p)
}

senders.iproc.design <- function(object, ...) {
    design <- object
    senders <- .Call("Riproc_design_senders", design)
    attributes(senders) <- attr(design, "attributes.senders")
    senders
}

receivers.iproc.design <- function(object, ...) {
    design <- object
    receivers <- .Call("Riproc_design_receivers", design)
    attributes(receivers) <- attr(design, "attributes.receivers")
    receivers
}

as.matrix.iproc.design <- function(object, sender, cursor = NULL, ...) {
    design <- object
    p <- ncol(design)
    mul(design, diag(p), sender, cursor)
}

mul.iproc.design <- function(x, y, sender, cursor = NULL, ...) {
    design <- x
    matrix <- as.matrix(y)
    storage.mode(matrix) <- "numeric"
    sender <- as.integer(sender)
    .Call("Riproc_design_mul", design, matrix, sender, cursor)
}

tmul.iproc.design <- function(x, y, sender, cursor = NULL, ...) {
    design <- x
    matrix <- as.matrix(y)
    storage.mode(matrix) <- "numeric"
    sender <- as.integer(sender)
    .Call("Riproc_design_tmul", design, matrix, sender, cursor)
}
