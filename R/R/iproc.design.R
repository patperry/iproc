
SPECIALS <- character(0)

var.interval <- function(name)
{
    specials <- get("SPECIALS", envir=parent.frame())
    assign("SPECIALS", c(specials, name), envir=parent.frame())

    function (interval) {
        interval <- as.numeric(interval, units="secs")
        if (length(interval) != 1L)
            stop("interval must be a scalar")
        if (is.na(interval))
            stop("interval is missing")
        if (interval <= 0)
            stop("interval is not positive")

        x <- list(interval = interval)
        class(x) <- c(name, "interval")
        x
    }
}

var.intervals <- function(name)
{
    specials <- get("SPECIALS", envir=parent.frame())
    assign("SPECIALS", c(specials, name), envir=parent.frame())

    function (...) {
        intervals <- c(...)
        intervals <- as.numeric(intervals, units="secs")
        if (length(intervals) == 0L)
            stop("must specify at least one interval")
        if (any(is.na(intervals)))
            stop("interval is missing")
        if (any(intervals <= 0))
            stop("interval is not positive")
        if (length(intervals) != length(unique(intervals)))
            stop("intervals must be unique")
        if (any(intervals != sort(intervals)))
            stop("intervals must be sorted in increasing order")

        x <- list(intervals = intervals)
        class(x) <- c(name, "intervals")
        x
    }
}

irecv <- var.interval("irecv")
isend <- var.interval("isend")
irecvtot <- var.interval("irecvtot")
isendtot <- var.interval("isendtot")

nrecv <- var.intervals("nrecv")
nsend <- var.intervals("nsend")
nrecvtot <- var.intervals("nrecvtot")
nsendtot <- var.intervals("nsendtot")


recv.frame <- function(formula, receiver.data = actor.data,
                       sender.data = actor.data,
                       actor.data = NULL, message.data = NULL,
                       sender.set = actor.set,
                       receiver.set = actor.set,
                       actor.set = NULL, contrasts = NULL, ...)
{
    call <- match.call()
    if (missing(message.data))
        message.data <- environment(formula)
    if (missing(actor.data))
        actor.data <- environment(formula)

    specials <- SPECIALS
    tm <- terms(formula, specials, data=receiver.data)
    v <- attr(tm, "variables")[-1L]
    nv <- nrow(attr(tm, "factors"))
    nt <- ncol(attr(tm, "factors"))

    # validate the formula
    if (any(attr(tm, "factors") == 2L))
        stop("interactions without main effects are not supported")
    if (attr(tm, "intercept") == 0L)
        warning("intercept term is irrelevant")
    if (attr(tm, "response") == 0L)
        stop("must specify a response")
    if (!is.null(attr(tm, "offset")) && attr(tm, "offset") != 0L)
        stop("offsets are not supported")


    # extract the response
    call.y <- v[[attr(tm, "response")]]
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
        stop("response must be a message object")


    # extract the model frame for the receivers
    mf.tvars <- list()
    is.trait <- rep(TRUE, nv)
    is.trait[[attr(tm, "response")]] <- FALSE

    term.labels <- attr(tm, "term.labels")
    special.ind <- unlist(attr(tm, "specials"))
    for (i in seq_along(term.labels)) {
        k <- match(i + 1L, special.ind)
        if (!is.na(k)) {
            call.i <- v[[i + 1L]]
            mf.tvars[[length(mf.tvars) + 1L]] <- eval(call.i, list(...))
            names(mf.tvars)[[length(mf.tvars)]] <- term.labels[[i]]
            is.trait[[i + 1L]] <- FALSE
        }
    }
    browser()
    call.x <- attr(tm, "variables")[c(1L, 1L + which(is.trait))]
    call.x[[1L]] <- quote(data.frame)
    data.x <- eval(call.x, receiver.data)
    mf.traits <- model.frame(~ . - 1, data.x)
    x.traits <- model.matrix(mf.traits, data.x, contrasts)
    # attr(x.traits, "assign") <- which(is.trait)[attr(x.traits, "assign")]

    browser()
}
