recv.model <- function(formula, message.data,
                       actor.data, receiver.data = actor.data,
                       sender.data = actor.data,
                       bipartite = FALSE, loops = FALSE,
                       actor.set = NULL, sender.set = NULL,
                       receiver.set = NULL,
                       method = "newton",
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

    y <- recv.frame.response(mf)
    if (is.null(y))
        stop("must specify a response")
    if (!inherits(y, "Mesg"))
        stop("response must be a message object")

    if (attr(mt, "intercept") == 0L)
        warning("intercept term is irrelevant")
    if (!is.null(attr(mt, "offset")) && attr(mt, "offset") != 0L)
        stop("offsets are not supported")

    if (nrow(attr(mt, "factors")) != ncol(attr(mt, "factors")) + 1L)
        stop("interactions are not implemented")


    time <- y$time
    sender <- y$sender
    receiver <- y$receiver
    attribute <- y$attribute

    sender.set <- if (is.null(mf$sender.set)) mf$actor.set
        else mf$sender.set
    receiver.set <- if (is.null(mf$receiver.set)) mf$actor.set
        else mf$receiver.set
    loops <- mf$loops

    factors <- attr(mt, "factors")
    types <- mf$types
    traits <- model.matrix(mf$traits, contrasts.arg=contrasts)
    specials <- mf$specials

    #model <- .Call("Riproc_recv_model", time, sender, receiver, attribute,
    #               sender.set, receiver.set, loops, factors, types, traits,
    #               specials)
    #model
}

