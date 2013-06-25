recv.model <- function(formula, message.data, receiver.data,
                       sender.data = receiver.data, n = nrow(receiver.data),
                       bipartite = FALSE, loops = FALSE, method = "newton",
                       contrasts = NULL, ...)
{
    call <- match.call()

    # compute the model frame
    mf <- match.call(expand.dots = TRUE)
    m <- match(c("contrasts", "method"), names(mf))
    m <- m[!is.na(m)]
    if (length(m) > 0)
        mf <- mf[-m]
    mf[[1L]] <- as.name("recv.frame")
    mf <- eval(mf, parent.frame())

    if (method == "recv.frame")
        return(mf)
    else if (method != "newton")
        warning(sprintf("method = '%s' is not supported. Using 'newton'",
                        method))

    y <- recv.response(mf)
    x <- recv.matrix(mf, contrasts=contrasts)

    if (is.null(y))
        stop("must specify a response")
    if (!inherits(y, "Mesg"))
        stop("response must be a message object")

    if (any(x$order > 1L))
        warning("interactions are not implemented")


    time <- y$time
    sender <- y$sender
    receiver <- y$receiver
    message.attr <- y$message.attr

    factors <- x$factors
    variables <- x$variables
    order <- x$order

    loops <- if (is.na(mf$loops)) TRUE else mf$loops

    nrecv <- n
    nsend <- if (bipartite) max(sender) else n

    browser()

    model <- .Call("Riproc_recv_model", time, sender, receiver, message.attr,
                   factors, variables, order, nsend, nrecv, loops)
    model
}

