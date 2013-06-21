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

    mt <- mf$terms
    factors <- get.factors(mt)

    y <- recv.frame.response(mf)
    if (is.null(y))
        stop("must specify a response")
    if (!inherits(y, "Mesg"))
        stop("response must be a message object")

    if (attr(mt, "intercept") == 0L)
        warning("intercept term is irrelevant")
    if (!is.null(attr(mt, "offset")) && attr(mt, "offset") != 0L)
        stop("offsets are not supported")

    if (nrow(factors) != ncol(factors) + 1L)
        stop("interactions are not implemented")


    time <- y$time
    sender <- y$sender
    receiver <- y$receiver
    message.attr <- y$message.attr

    loops <- if (is.na(mf$loops)) TRUE else mf$loops

    types <- mf$types
    if (!is.null(mf$traits)) {
        traits <- model.matrix(mf$traits, contrasts.arg=contrasts)
    } else {
        traits <- data.frame(rep(1, n))
        names(traits) <- "(Intercept)"
        attr(traits, "assign") <- 0L
    }
    specials <- mf$specials

    nrecv <- n
    nsend <- if (bipartite) max(sender) else n

    browser()

    model <- .Call("Riproc_recv_model", time, sender, receiver, message.attr,
                   nsend, nrecv, loops, factors, types, traits,
                   specials)
    model
}

