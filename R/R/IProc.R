IProc <- function(time, sender, receiver, attribute = NULL, data = NULL)
{
    if (missing(time))
        stop("Must have a time argument")
    if (missing(sender))
        stop("Must have a sender argument")
    if (missing(receiver))
        stop("Must have a receiver argument")

    time <- eval(call("with", substitute(data), substitute(time)))
    sender <- eval(call("with", substitute(data), substitute(sender)))
    receiver <- eval(call("with", substitute(data), substitute(receiver)))
    attribute <- eval(call("with", substitute(data), substitute(attribute)))


    # validate and normalize time
    if (!is.numeric(unclass(time)) && !is.integer(time))
        stop("Time variable is not numeric, integer, or POSIXct")
    if (any(is.na(time)))
        stop("Time variable contains missing values")
    if (any(is.nan(time)))
        stop("Time variable contains NaN values")

    storage.mode(time) <- "numeric"


    # validate and normalize sender
    if (!all(sender == as.integer(sender))) # fails on NaN, NA, Inf
        stop("Sender variable contains non-integer values")
    if (length(sender) != length(time))
        stop(sprintf("Number of senders is %d, should equal %d (number of times)",
             length(sender), length(time)))

    sender <- as.integer(sender)
    sender.orig <- sort.int(unique(sender))
    sender <- match(sender, sender.orig)


    # validate and normalize receiver
    receiver <- as.list(receiver)
    if (!all(unlist(lapply(receiver, as.integer)) == unlist(receiver)))
        stop("Receiver variable contains non-integer values")

    if (length(receiver) != length(time))
        stop(sprintf("Number of receivers is %d, should equal %d (number of times)",
             length(receiver), length(time)))

    receiver <- lapply(receiver, as.integer)
    receiver.orig <- sort.int(unique(unlist(receiver)))
    receiver <- lapply(receiver, match, table=receiver.orig)


    # validate and normalize attribute
    attribute <- as.vector(attribute)

    if (!is.null(attribute)) {
        if (length(attribute) != length(time))
            stop(sprintf("Number of attributes is %d, should equal %d (number of times)",
                         length(attribute), length(time)))
    }


    # sort observations by time
    id <- order(time)
    time <- time[id]
    sender <- sender[id]
    receiver <- receiver[id]
    attribute <- attribute[id]

    z <- data.frame(time = I(time), sender = sender, receiver = I(receiver))
    if (!is.null(attribute))
        z$attribute <- I(attribute)

    attr(z, "sender.orig") <- sender.orig
    attr(z, "receiver.orig") <- receiver.orig
    class(z) <- c("IProc", class(z))

    z
}


"[.IProc" <- function(x, i, j, drop=FALSE) {
    if (missing(j)) {
        id <- attr(x, "id")
        id.orig <- attr(x, "id.orig")
        time.attr <- attr(x, "time")
        sender.orig <- attr(x, "sender.orig")
        receiver.orig <- attr(x, "receiver.orig")
        ctemp <- class(x)
        class(x) <- "data.frame"
    }
}
