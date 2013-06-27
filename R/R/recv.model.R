recv.model <- function(formula, message.data, receiver.data,
                       sender.data = receiver.data, n = nrow(receiver.data),
                       bipartite = FALSE, loops = FALSE, method = "newton",
                       contrasts = NULL, skip = recv.horizon, ...)
{
    call <- match.call()

    # compute the model frame
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("contrasts", "method", "skip"), names(mf))
    m <- m[!is.na(m)]
    if (length(m) > 0)
        mf <- mf[-m]
    mf[[1L]] <- as.name("recv.frame")
    mf <- eval(mf, parent.frame())
    recv.horizon <- recv.horizon(mf)

    if (method == "recv.frame")
        return(mf)
    else if (method != "newton")
        warning(sprintf("method = '%s' is not supported. Using 'newton'",
                        method))

    if (!(skip >= 0))
        stop("'skip' must be non-negative")

    y <- recv.response(mf)
    x <- recv.matrix(mf, contrasts=contrasts)

    if (is.null(y))
        stop("must specify a response")
    if (!inherits(y, "Mesg"))
        stop("response must be a message object")


    time <- y$time
    sender <- y$sender
    receiver <- y$receiver

    factors <- x$factors
    variables <- x$variables

    loops <- if (is.na(mf$loops)) TRUE else mf$loops

    nrecv <- n
    nsend <- if (bipartite) max(sender) else n

    model <- .Call("Riproc_recv_model", time, sender, receiver,
                   factors, variables, nsend, nrecv, loops, skip)
    model
}

