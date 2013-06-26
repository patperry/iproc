
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
        attr(x, "horizon") <- interval
	attr(x, "labels") <- paste(name, "(", interval, ")", sep="")
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
        attr(x, "horizon") <- max(intervals)
	attr(x, "labels") <- paste(name, "(", intervals, ")", sep="")
        x
    }
}

var.intervals2 <- function(name, type, bipartite = FALSE)
{
    specials <- get("SPECIALS", envir=parent.frame())
    assign("SPECIALS", c(specials, name), envir=parent.frame())

    if (bipartite) {
        bspecials <- get("BSPECIALS", envir=parent.frame())
        assign("BSPECIALS", c(bspecials, name), envir=parent.frame())
    }

    function (intervals1, intervals2) {
        intervals1 <- as.numeric(intervals1, units="secs")
        intervals2 <- as.numeric(intervals2, units="secs")
        for (intervals in list(intervals1, intervals2)) {
            if (length(intervals) == 0L)
                stop("must specify at least one interval in each argument")
            if (any(is.na(intervals)))
                stop("interval is missing from")
            if (any(intervals <= 0))
                stop("interval is not positive")
            if (length(intervals) != length(unique(intervals)))
                stop("intervals must be unique")
            if (any(intervals != sort(intervals)))
                stop("intervals must be sorted in increasing order")
        }

        x <- list(intervals1 = intervals1, intervals2 = intervals2)
        class(x) <- c(name, type, "intervals2")
        attr(x, "horizon") <- max(intervals1, intervals2)
	attr(x, "labels") <- paste(name, "(",
				   rep(intervals1, times=length(intervals2)),
				   ",",
				   rep(intervals2, each=length(intervals1)),
				   ")", sep="")
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

nrecv2 <- var.intervals2("nrecv2", "dyad")
nsend2 <- var.intervals2("nsend2", "dyad")
nsib <- var.intervals2("nsib", "dyad")
ncosib <- var.intervals2("ncosib", "dyad")



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
            str <- paste(deparse(ret), collapse=" ")
            ret <- as.formula(paste("~ 1 +", str))[[2L]]
        }
    } else {
        ret <- object
        args <- lapply(object[-1L], expand.s, specials=specials, data=data)
        ret[-1L] <- args
    }

    environment(ret) <- environment(object)
    ret
}


recv.frame.variable <- function(object, specials, data, sender.data,
                                envir=parent.frame(), ...)
{
    if (is.call(object) && object[[1L]] == "s") {
        if (length(object) != 2L)
            stop("Too many arguments to 's'")

        var <- recv.frame.variable(object[[2L]], specials, sender.data, NULL,
                                   envir, ...)
        name <- paste(deparse(object[[2L]]), collapse=" ")

        if (inherits(var, "dyad"))
            stop(sprintf("term 's(%s)' is invalid; '%s' is a dyad-specific variable",
                         name, name))

        class(var) <- c("send", class(var))

    } else if (is.call(object)
               && as.character(object[[1L]]) %in% specials) {
        var <- eval(object, list(...))
        class(var) <- c(class(var), "special")
    } else {
        data <- eval(call("data.frame", object), data)

        env <- as.environment(data)
        parent.env(env) <- envir

        mf <- model.frame(~ ., data, na.action=na.fail)
        environment(mf) <- env

        var <- mf
        #class(var) <- c("trait", class(var))
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
    response <- attr(tm, "response")

    # validate the formula
    special.names <- names(attr(tm, "specials"))[!sapply(attr(tm, "specials"),
                                                         is.null)]
    if (any(attr(tm, "factors") == 2L))
        stop("interactions without main effects are not supported")
    if (any(attr(tm, "factors")[response,-response] == 1L))
        stop("interactions between the response and the predictors are not supported")

    for (s in setdiff(specials, bspecials)) {
        if (bipartite && s %in% special.names)
            stop(sprintf("the special term '%s' is invalid for bipartite processes", s))
    }

    v <- attr(tm, "variables")[-1L]
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

            variables[[i]] <- y
        } else {
            rv <- recv.frame.variable(v[[i]], specials, receiver.data,
                                      sender.data, ...)
            name <- attr(rv, "name")
            attr(rv, "name") <- NULL
            variables[[i]] <- rv
        }
    }

    send <- sapply(variables, function(x) inherits(x, "send"))
    factors <- attr(tm, "factors")
    term.labels <- attr(tm, "term.labels")
    order <- attr(tm, "order")

    all.send <- colSums(factors * !send) == 0
    factors <- factors[,!all.send,drop=FALSE]
    term.labels <- term.labels[!all.send]
    order <- order[!all.send]

    frame <- list(formula = formula(tm), factors = factors,
                  term.labels = term.labels, variables = variables,
                  order = order, response = response,
                  n = n, bipartite = bipartite, loops = loops)
    class(frame) <- "recv.frame"
    frame
}


recv.response <- function(frame)
{
    if (!inherits(frame, "recv.frame"))
        stop("invalid 'frame' argument")

    response <- frame$response
    y <- frame$variables[[response]]
    y
}


recv.weights <- function(frame)
{
    stop("not implemented")
}


recv.horizon <- function(frame)
{
    if (!inherits(frame, "recv.frame"))
        stop("invalid 'frame' argument")

    hs <- lapply(frame$variables, function(x) attr(x, "horizon", TRUE))
    horizon <- max(c(0, unlist(hs)))
    horizon
}




