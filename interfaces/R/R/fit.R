

iproc.fit <- function(model0, messages, penalty, control = list()) {
    penalty <- as.numeric(penalty)
    control <- do.call("iproc.control", control)
    model <- .Call("Riproc_fit", model0, messages, penalty, control$reltol,
                   control$abstol, control$maxit, control$trace,
                   control$report)
    model
}

iproc.control <- function(reltol = .Machine$double.eps,
                          abstol = sqrt(reltol),
                          maxit = 1000L,
                          trace = FALSE,
                          report = 10L)
{
    if (!is.numeric(reltol) || reltol <= 0)
        stop("value of 'reltol' must be > 0")
    if (!is.numeric(abstol) || abstol <= 0)
        stop("value of 'abstol' must be > 0")
    if (!(is.integer(maxit) || is.numeric(maxit)) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    if (!(is.integer(report) || is.numeric(report)) || report <= 0)
        stop("value of 'report' must be > 0")

    trace <- as.logical(trace)
    if (is.na(trace))
        trace <- FALSE
    
    list(reltol = reltol, abstol = abstol, maxit = as.integer(maxit),
         trace = trace, report = as.integer(report))
}
