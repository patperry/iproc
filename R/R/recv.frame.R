
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


recv.frame <- function(formula, message.data = NULL, actor.data = NULL,
                       sender.data = actor.data, receiver.data = actor.data,
                       bipartite = is.null(actor.set), loops = FALSE,
                       actor.set = NULL, sender.set = NULL,
                       receiver.set = NULL, ...)
{
    call <- match.call()
    if (missing(message.data))
        message.data <- environment(formula)
    if (missing(actor.data))
        actor.data <- environment(formula)

    if (!(is.logical(bipartite) && length(bipartite) == 1L))
        stop("'bipartite' must be a logical scalar")
    if (!(is.logical(loops) && length(loops) == 1L))
        stop("'loops' must be a logical scalar")


    specials <- SPECIALS
    tm <- terms(formula, specials, data=receiver.data)
    v <- attr(tm, "variables")[-1L]
    nv <- nrow(attr(tm, "factors"))
    nt <- ncol(attr(tm, "factors"))
    special.ind <- unlist(attr(tm, "specials"))
    special.names <- names(attr(tm, "specials"))[!sapply(attr(tm, "specials"), is.null)]


    # validate the formula
    if (any(attr(tm, "factors") == 2L))
        stop("interactions without main effects are not supported")
    for (s in c("irecvtot", "nrecvtot", "irecv", "nrecv",
                "nsend2", "nrecv2", "nsib", "ncosib")) {
        if (bipartite && s %in% special.names)
            stop(sprintf("the special term '%s' is invalid for bipartite processes", s))
    }


    # extract the response (if one exists)
    if (attr(tm, "response") != 0L) {
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

        if (is.Mesg(y)) {
            if (bipartite) {
                sender.set <- attr(y, "sender.set")
                receiver.set <- attr(y, "receiver.set")
            } else {
                actor.set <- sort(unique(c(attr(y, "sender.set"),
                                           attr(y, "receiver.set"))))
            }
        }
    } else {
        y <- NULL
    }

    if (bipartite) {
        if (is.null(sender.set))
            stop("must specify sender.set")
        if (length(sender.set) == 0L)
            stop("must have at least one sender")
        if (is.null(receiver.set))
            stop("must specify receiver.set")
        if (length(receiver.set) == 0L)
            stop("must have at least one receiver")
        if (!is.null(actor.set))
            stop("cannot specify actor.set for bipartite processes")
    } else {
        if (is.null(actor.set))
            stop("must specify actor.set")
        if (length(actor.set) == 0L)
            stop("must have at least one actor")
        if (!loops && length(actor.set) == 1L)
            stop("must have at least two actors")
        if (!is.null(sender.set))
            stop("cannot specify sender.set for bipartite processes")
        if (!is.null(receiver.set))
            stop("cannot specify receiver.set for bipartite processes")

    }


    # extract the model frame for the specials and traits
    mf.specials <- list()

    names <- rownames(attr(tm, "factors"))
    is.trait <- rep(TRUE, nv)
    if (attr(tm, "response") != 0L) {
        i <- attr(tm, "response")
        is.trait[[i]] <- FALSE
    }

    for (i in seq_along(names)) {
        k <- match(i, special.ind)
        if (!is.na(k)) {
            call.i <- v[[i]]
            mf.specials[[length(mf.specials) + 1L]] <- eval(call.i, list(...))
            names(mf.specials)[[length(mf.specials)]] <- names[[i]]
            is.trait[[i]] <- FALSE
        }
    }

    call.x <- attr(tm, "variables")[c(1L, 1L + which(is.trait))]
    call.x[[1L]] <- quote(data.frame)
    data.x <- eval(call.x, receiver.data)
    mf.traits <- model.frame(~ . - 1, data.x)

    frame <- list(formula = formula(tm),
                  terms = tm, response = y, traits = mf.traits,
                  traits.data = data.x,
                  specials = mf.specials, sender.set = sender.set,
                  receiver.set = receiver.set, actor.set = actor.set,
                  bipartite = bipartite, loops = loops)
    class(frame) <- "recv.frame"
    frame
}


recv.frame.response <- function(frame)
{
    if (!inherits(frame, "recv.frame"))
        stop("invalid 'frame' argument")
    frame$response
}

recv.model.weights <- function(frame)
{
    stop("not implemented")
}


model.matrix.recv.frame <- function(object, time, sender, ...)
{
    stop("not implemented")
}
