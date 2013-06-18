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


    # validate and normalize receiver
    receiver <- as.list(receiver)
    if (!all(as.integer(unlist(receiver)) == unlist(receiver)))
        stop("Receiver variable contains non-integer values")
    if (length(receiver) != length(time))
        stop(sprintf("Number of receivers is %d, should equal %d (number of times)",
             length(receiver), length(time)))

    receiver <- lapply(receiver, as.integer)


    # validate and normalize attribute
    if (!is.null(attribute)) {
        if (length(attribute) != length(time))
            stop(sprintf("Number of attributes is %d, should equal %d (number of times)",
                         length(attribute), length(time)))
    }

    attribute <- as.vector(attribute)


    # convert ids to compact form
    sender.orig <- sort.int(unique(sender))
    sender <- match(sender, sender.orig)
    receiver.orig <- sort.int(unique(unlist(receiver)))
    receiver <- lapply(receiver, match, table=receiver.orig)

    # sort observations by time
    id <- order(time)
    id.orig <- rank(unclass(time)) # 'unclass' is a work-around for a bug in R-2.15
    time <- time[id]
    sender <- sender[id]
    receiver <- receiver[id]
    attribute <- attribute[id]

    z <- data.frame(time = I(time), sender = sender, receiver = I(receiver))
    if (!is.null(attribute))
        z$attribute <- I(attribute)

    attr(z, "id") <- id
    attr(z, "id.orig") <- id.orig
    attr(z, "sender.orig") <- sender.orig
    attr(z, "receiver.orig") <- receiver.orig
    class(z) <- c("IProc", class(z))

    z
}


print.IProc <- function(x, ...)
{
    str <- as.character.IProc(x)
    n <- length(str)

    if (n > 0) {
        i <- seq_len(n)
        d <- 1L + floor(log10(i))
        dmax <- d[n]
        spaces <- "                               " # 31 spaces
        pad <- sapply(dmax - d, substr, x=spaces, start=1)

        cat(paste(pad, "[", i, "] ", str, sep=""), sep="\n")

    } else {
        cat("IProc(0)")
    }

    invisible(str)
}


as.character.IProc <- function(x, ...)
{
    tw <- if(inherits(x$time, "POSIXct")) 80 else NULL
    nr <- sapply(x$receiver, length)
    paste(format(x$time, width=tw), " : ",
          x$sender, " -> ",
          ifelse(nr > 1, "[ ", ""),
            format(x$receiver),
          ifelse(nr > 1, " ]", ""),
          sep="")
}


"[.IProc" <- function(x, i, j, drop = FALSE)
{
    id <- attr(x, "id")
    id.orig <- attr(x, "id.orig")
    sender.orig <- attr(x, "sender.orig")
    receiver.orig <- attr(x, "receiver.orig")

    i <- if (missing(i))
        seq_len(nrow(x))
    else id.orig[i]

    # If only 1 subscript is given, the result will still be an IProc object,
    # and the drop argument is ignored.
    if (missing(j)) {
        ctemp <- class(x)
        class(x) <- "data.frame"
        x <- x[i,, drop=FALSE]
        class(x) <- ctemp
        attr(x, "id") <- id[i]
        attr(x, "id.orig") <- i
        attr(x, "sender.orig") <- sender.orig
        attr(x, "receiver.orig") <- receiver.orig
        x
    }
    else { # return  a matrix or vector
        class(x) <- "data.frame"
        x$sender <- sender.orig[x$sender]
        x$receiver <- lapply(x$receiver, function(r) receiver.orig[r])
        NextMethod("[")
    }
}

"$.IProc" <- function(x, name)
{
    x[,name]
}

"$<-.IProc" <- function(x, name, value) stop("Cannot reassign intersection process values")


"[[.IProc" <- function(x, i, j, ..., exact = TRUE)
{
    x <- as.matrix(x)
    NextMethod("[[")
}

"[[<-.IProc" <- function(x, i, j, ..., value) stop("Cannot reassign interaction process values")


is.na.IProc <- function(x) rep(FALSE, nrow(x))

Math.IProc <- function(...) stop("Invalid operation on an interaction process")
Ops.IProc <- function(...) stop("Invalid operation on an interaction process")
Summary.IProc <- function(...) stop("Invalid operation on an interaction process")
is.IProc <- function(x) inherits(x, "IProc")

as.matrix.IProc <- function(x, ...) {
    x[seq_len(nrow(x)), seq_len(ncol(x)), drop=FALSE]
}


summary.IProc <- function(object, ...) {
    ans <- list()
    ans$n <- nrow(object)
    ans$nsender <- length(attr(object, "sender.orig"))
    ans$nreceiver <- length(attr(object, "receiver.orig"))
    class(ans) <- "summary.IProc"

    ans
}


print.summary.IProc <- function(x, ...) {
    cat(sprintf("%d messages (%d senders, %d receivers)",
                x$n, x$nsender, x$nreceiver))
    cat("\n")
    invisible(x)
}
