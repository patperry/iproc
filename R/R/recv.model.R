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
    perm <- model$perm
    names <- model$names

    model$coefficients <- model$coefficients[perm]
    names(model$coefficients) <- names

    model$constraints <- model$constraints[,perm,drop=FALSE]
    colnames(model$constraints) <- names

    model$score <- model$score[perm]
    names(model$score) <- names

    model$imat <- model$imat[perm,perm,drop=FALSE]
    rownames(model$imat) <- names
    colnames(model$imat) <- names

    model$perm <- NULL
    model$names <- NULL


    model$call <- call
    class(model) <- "recv.model"
    model
}


print.recv.model <- function(x, digits = max(3, getOption("digits") - 3), ...)
{
    print.glm(x, digits, ...)
}


vcov.recv.model <- function(object, ...)
{
    imat <- object$imat
    constr <- object$constraints

    n <- nrow(imat)
    k <- nrow(constr)
    N <- n + k

    i <- seq_len(n)
    j <- n + seq_len(k)

    K <- matrix(NA, N, N)
    K[i, i] <- imat
    K[i, j] <- t(constr)
    K[j, i] <- constr
    K[j, j] <- 0

    cov <- solve(K)
    cov <- cov[i, i, drop=FALSE]
    rownames(cov) <- rownames(imat)
    colnames(cov) <- colnames(imat)
    cov
}
