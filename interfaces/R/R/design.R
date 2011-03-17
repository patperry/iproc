
iproc.design.default <- function(senders, receivers = senders, recip.intervals = NULL,
                                 data.senders, data.receivers = data.senders, ...)
{
    if (!inherits(senders, "actors")) {
        if (missing(data.senders)) {
            senders <- actors(senders)
        } else {
            senders <- actors(senders, data.senders)
        }
    }
    
    if (!inherits(receivers, "actors")) {
        if (missing(data.receivers)) {
            receivers <- actors(receivers)
        } else {
            receivers <- actors(receivers, data.receivers)
        }
    }

    recip.intervals.orig <- as.numeric(recip.intervals)
    if (any(is.na(recip.intervals.orig) | (recip.intervals.orig <= 0)))
        stop("'recip.intervals' contains a nonpositive or NA value")
    
    recip.intervals <- unique(sort(recip.intervals.orig))
    if (!identical(recip.intervals, recip.intervals.orig))
        stop("'recip.intervals' not unique and sorted in ascending order")

    design <- .Call("Riproc_design_new", senders, receivers, recip.intervals)
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
