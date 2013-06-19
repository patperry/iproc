
irecv <- function(interval)
{
    interval <- as.numeric(interval, units="secs")
    if (length(interval) != 1L)
        stop("Interval must be a scalar")
    if (is.na(interval))
        stop("Interval is missing")
    if (interval <= 0)
        stop("Interval is not positive")

    list(irecv=interval)
}


iproc.design <- function(formula, receiver.data = actor.data,
                         sender.data = actor.data,
                         actor.data = NULL, message.data = NULL,
                         sender.set = actor.set,
                         receiver.set = actor.set,
                         actor.set = NULL, ...)
{
    call <- match.call()
    if (missing(message.data))
        message.data <- environment(formula)
    if (missing(actor.data))
        actor.data <- environment(formula)

    tm <- terms(formula, data=message.data, specials="irecv")
    v <- attr(tm, "variables")

    if (attr(tm, "offset") != 0L)
        top("Offsets are not supported")

    vx <- v[[-(1L + attr(tm, "response"))]]

    # construct the response object
    if (attr(tm, "response") == 0L)
        stop("Must specify a response")
    call.y <- v[[1L + attr(tm, "response")]]
    ny <- length(call.y)

    if (!is.null(sender.set)) {
        ny <- ny + 1L
        call.y[[ny]] <- substitute(sender.set)
        names(call.y)[[ny]] <- "sender.set"
    }
    if (!is.null(receiver.set)) {
        ny <- ny + 1L
        call.y[[ny]] <- substitute(receiver.set)
        names(call.y)[[ny]] <- "receiver.set"
    }
    y <- eval(call.y, message.data)
    if (!inherits(y, "Mesg"))
        stop("Response must be a message object")

    browser()
}
