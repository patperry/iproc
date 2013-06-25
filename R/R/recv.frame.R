
SPECIALS <- character(0)
BSPECIALS <- character(0)

var.interval <- function(name, type, bipartite = FALSE)
{
    specials <- get("SPECIALS", envir=parent.frame())
    assign("SPECIALS", c(specials, name), envir=parent.frame())

    if (bipartite) {
        bspecials <- get("BSPECIALS", envir=parent.frame())
        assign("BSPECIALS", c(bspecials, name), envir=parent.frame())
    }

    function (interval) {
        interval <- as.numeric(interval, units="secs")
        if (length(interval) != 1L)
            stop("interval must be a scalar")
        if (is.na(interval))
            stop("interval is missing")
        if (interval <= 0)
            stop("interval is not positive")

        x <- list(interval = interval)
        class(x) <- c(name, type, "interval")
        x
    }
}

var.intervals <- function(name, type, bipartite = FALSE)
{
    specials <- get("SPECIALS", envir=parent.frame())
    assign("SPECIALS", c(specials, name), envir=parent.frame())

    if (bipartite) {
        bspecials <- get("BSPECIALS", envir=parent.frame())
        assign("BSPECIALS", c(bspecials, name), envir=parent.frame())
    }

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
        class(x) <- c(name, type, "intervals")
        x
    }
}

irecv <- var.interval("irecv", "dyad", bipartite=TRUE)
nrecv <- var.intervals("nrecv", "dyad", bipartite=TRUE)
irecvtot <- var.interval("irecvtot", "actor", bipartite=TRUE)
nrecvtot <- var.intervals("nrecvtot", "actor", bipartite=TRUE)


isend <- var.interval("isend", "dyad")
isendtot <- var.interval("isendtot", "actor")
nsend <- var.intervals("nsend", "dyad")
nsendtot <- var.intervals("nsendtot", "actor")



add.s <- function(object, specials = NULL)
{
    if (is.name(object)
        || (is.call(object) && as.character(object[[1L]]) %in% specials)) {
        ret <- call("s", object)
    } else {
        ret <- object
        args <- lapply(object[-1L], add.s, specials=specials)
        ret[-1] <- args
    }
    ret
}


expand.s <- function(object, specials = NULL, data = NULL) {
    if (missing(data))
        data <- environment(object)

    if (!is.call(object)) {
        ret <- object
    } else if (object[[1L]] == quote(s)) {
        fmla <- object
        fmla[[1L]] <- quote(`~`)
        tm <- eval(call("terms", fmla, specials=specials, data=data))
        labels <- attr(tm, "term.labels")
        parse <- as.formula(paste("~", paste(labels, collapse=" + ")))[[2]]
        ret <- add.s(parse, specials=specials)
        if (attr(tm, "intercept") == 1L) {
            ret <- as.formula(paste("~ 1 +", deparse(ret)))[[2L]]
        }
    } else {
        ret <- object
        args <- lapply(object[-1L], expand.s, specials=specials, data=data)
        ret[-1L] <- args
    }

    ret
}


recv.variable <- function(object, specials, data, sender.data, ...)
{
    if (is.call(object) && object[[1L]] == "s") {
        if (length(object) != 2L)
            stop("Too many arguments to 's'")

        var <- recv.variable(object[[2L]], specials, sender.data, NULL, ...)
        name <- deparse(object[[2L]])

        if (inherits(var, "dyad"))
            stop(sprintf("term 's(%s)' is invalid; '%s' is a dyad-specific variable",
                         name, name))

        if (inherits(var, "trait"))
            colnames(var) <- paste("s(", colnames(var), ")", sep="")

        class(var) <- c("send", class(var))
    } else if (is.call(object)
               && as.character(object[[1L]]) %in% specials) {
        var <- eval(object, list(...))
        class(var) <- c("special", class(var))
    } else {
        data <- eval(call("data.frame", object), data)
        mm <- model.matrix(~ ., data)
        int <- match("(Intercept)", colnames(mm))
        if (!is.na(int)) {
            assign <- attr(mm, "assign")
            contrasts <- attr(mm, "contrasts")
            var <- mm[,-int,drop=FALSE]
            attr(var, "contrasts") <- contrasts[[1L]]
            class(var) <- c("trait", class(mm))
        }
    }

    var
}



recv.frame <- function(formula, message.data = NULL, receiver.data = NULL,
                       sender.data = receiver.data, n = nrow(receiver.data),
                       bipartite = FALSE, loops = FALSE, ...)
{
    if (missing(message.data))
        message.data <- environment(formula)
    if (missing(receiver.data))
        receiver.data <- environment(formula)

    if (missing(n) && missing(receiver.data))
        stop("must specify 'n'")
    if (!(length(n) == 1L && n == as.integer(n) && n > 0L))
        stop("'n' must be a positive integer scalar")
    if (!(is.logical(bipartite) && length(bipartite) == 1L))
        stop("'bipartite' must be a logical scalar")
    if (!(is.logical(loops) && length(loops) == 1L))
        stop("'loops' must be a logical scalar")
    if (!bipartite && !loops && n == 1L)
        stop("must have at least two receivers")
    if (bipartite && !missing(loops) && !is.null(loops) && !is.na(loops)) {
        warning("'loops' argument is ignored for bipartite networks")
        loops <- NA
    }

    specials <- SPECIALS
    bspecials <- BSPECIALS
    formula <- expand.s(formula, specials, data=sender.data)
    tm <- terms(formula, specials, data=receiver.data)
    v <- attr(tm, "variables")[-1L]

    response <- attr(tm, "response")
    variables <- list()

    for (i in seq_along(v)) {
        if (i == response) {
            y <- eval(v[[i]], message.data)

            if (is.Mesg(y)) {
                if (!bipartite && max(y$sender) > n)
                    stop("message sender is out of range")
                if (max(unlist(y$receiver)) > n)
                    stop("message receiver is out of range")
            }

            class(y) <- c("response", class(y))
            variables[[i]] <- y
        } else {
            rv <- recv.variable(v[[i]], specials, receiver.data, sender.data, ...)
            name <- attr(rv, "name")
            attr(rv, "name") <- NULL
            variables[[i]] <- rv
        }
        names(variables)[[i]] <- deparse(v[[i]])
    }

    browser()

    factors <- get.factors(tm)
    nv <- nrow(factors)
    nt <- ncol(factors)
    special.ind <- unlist(attr(tm, "specials"))
    special.names <- names(attr(tm, "specials"))[!sapply(attr(tm, "specials"), is.null)]


    # validate the formula
    if (any(factors == 2L))
        stop("interactions without main effects are not supported")

    for (s in setdiff(specials, bspecials)) {
        if (bipartite && s %in% special.names)
            stop(sprintf("the special term '%s' is invalid for bipartite processes", s))
    }


    # extract the response (if one exists)
    if (attr(tm, "response") != 0L) {
        call.y <- v[[attr(tm, "response")]]
        y <- eval(call.y, message.data)

    } else {
        y <- NULL
    }



    # extract the model frame for the specials and traits
    mf.specials <- list()

    names <- rownames(factors)
    types <- rep(NA, nv)

    for (i in seq_along(names)) {
        k <- match(i, special.ind)
        if (!is.na(k)) {
            call.i <- v[[i]]
            mf.specials[[length(mf.specials) + 1L]] <- eval(call.i, list(...))
            names(mf.specials)[[length(mf.specials)]] <- names[[i]]
            types[[i]] <- "special"
        } else if (i == attr(tm, "response")) {
            types[[i]] <- "response"
        } else {
            types[[i]] <- "trait"
        }
    }

    if (any(types == "trait")) {
        fmla.x <- ~ .
        call.x <- attr(tm, "variables")[c(1L, 1L + which(types == "trait"))]
        call.x[[1L]] <- quote(data.frame)
        data.x <- eval(call.x, receiver.data)

        env.x <- as.environment(data.x)
        parent.env(env.x) <- parent.frame()

        mf.traits <- model.frame(fmla.x, data.x)
        environment(mf.traits) <- env.x
    } else {
        mf.traits <- NULL
    }


    frame <- list(formula = formula(tm), terms = tm, types = types,
                  response = y, traits = mf.traits, specials = mf.specials,
                  n = n, bipartite = bipartite, loops = loops)
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
