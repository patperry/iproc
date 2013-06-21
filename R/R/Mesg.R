Mesg <- function(time, sender, receiver, data, message.attr = NULL)
{
    if (missing(time))
        stop("must have a time argument")
    if (missing(sender))
        stop("must have a sender argument")
    if (missing(receiver))
        stop("must have a receiver argument")
    if (missing(data))
        data <- parent.frame()

    time <- eval(substitute(time), data)
    sender <- eval(substitute(sender), data)
    receiver <- eval(substitute(receiver), data)
    if (!missing(message.attr))
        message.attr <- eval(substitute(message.attr), data)


    # validate and normalize time
    if (!is.numeric(unclass(time)) && !is.integer(time))
        stop("time variable is not numeric, integer, or POSIXct")
    if (any(is.na(time)))
        stop("time variable contains missing values")
    if (any(is.nan(time)))
        stop("time variable contains NaN values")
    if (any(is.infinite(time)))
        stop("time variable contains infinite values")

    storage.mode(time) <- "numeric"


    # validate and normalize sender
    if (!all(sender == as.integer(sender))) # fails on NaN, NA, Inf
        stop("sender variable contains non-integer values")
    if (!all(sender > 0L))
        stop("sender variable contains non-positive values")
    if (length(sender) != length(time))
        stop(sprintf("number of senders is %d, should equal %d (number of times)",
             length(sender), length(time)))

    sender <- as.integer(sender)


    # validate and normalize receiver
    receiver <- as.list(receiver)
    if (!all(as.integer(unlist(receiver)) == unlist(receiver)))
        stop("receiver variable contains non-integer values")
    if (!all(unlist(receiver) > 0L))
        stop("receiver variable contains non-positive values")
    if (length(receiver) != length(time))
        stop(sprintf("number of receivers is %d, should equal %d (number of times)",
             length(receiver), length(time)))
    if (length(unlist(lapply(receiver, unique))) != length(unlist(receiver)))
        stop("receiver variable contains repeated entries (on same message)")

    receiver <- lapply(receiver, as.integer)


    # validate and normalize message.attr
    if (!is.null(message.attr)) {
        if (length(message.attr) != length(time))
            stop(sprintf("number of attributes is %d, should equal %d (number of times)",
                         length(message.attr), length(time)))
    }


    if (is.unsorted(time))
        warning("time is not sorted in ascending order; efficiency may be suboptimal")

    z <- data.frame(time = time, sender = sender)
    z$receiver <- receiver # need a separate assignment in case receiver is a list
    if (!is.null(message.attr))
        z$message.attr <- message.attr

    class(z) <- c("Mesg", class(z))

    z
}


print.Mesg <- function(x, ...)
{
    str <- as.character.Mesg(x)
    n <- length(str)
    names <- rownames(x)

    if (is.null(names)) {
        if (n > 0) {
            i <- seq_len(n)
            d <- 1L + floor(log10(i))
            dmax <- d[n]
            spaces <- "                               " # 31 spaces
            pad <- sapply(dmax - d, substr, x=spaces, start=1)

            cat(paste(pad, "[", i, "] ", str, sep=""), sep="\n")

        } else {
            cat("Mesg(0)")
        }
    } else {
        cat(paste(format(names), str, sep=" "), sep="\n")
    }

    invisible(str)
}


as.character.Mesg <- function(x, ...)
{
    tw <- if(inherits(x$time, "POSIXct")) 80 else NULL
    nr <- sapply(x$receiver, length)
    str <- paste(format(x$time, width=tw), " : ",
                 x$sender, " -> ",
                 ifelse(nr > 1, "[ ", ""),
                 format(x$receiver),
                 ifelse(nr > 1, " ]", ""),
                 sep="")
    names(str) <- rownames(x)
    str
}


is.na.Mesg <- function(x)
{
    res <- data.frame(time = is.na(x$time),
                      sender = is.na(x$sender),
                      receiver = sapply(x$receiver, function(r)
                                        is.null(r) || length(r) == 0L || any(is.na(r))))
    if (!is.null(x$message.attr)) {
        res$message.attr <- sapply(x$message.attr, is.na)
    }

    res
}

Math.Mesg <- function(...) stop("invalid operation on a message")
Ops.Mesg <- function(...) stop("invalid operation on a message")
Summary.Mesg <- function(...) stop("invalid operation on a message")
is.Mesg <- function(x) inherits(x, "Mesg")

as.matrix.Mesg <- function(x, ...) {
    x[seq_len(nrow(x)), seq_len(ncol(x)), drop=FALSE]
}


summary.Mesg <- function(object, ...) {
    ans <- list()
    ans$n <- nrow(object)
    ans$nsender <- max(object$sender)
    ans$nreceiver <- max(unlist(object$receiver))
    class(ans) <- "summary.Mesg"

    ans
}


print.summary.Mesg <- function(x, ...) {
    cat(sprintf("%d messages (%d senders, %d receivers)",
                x$n, x$nsender, x$nreceiver))
    cat("\n")
    invisible(x)
}

